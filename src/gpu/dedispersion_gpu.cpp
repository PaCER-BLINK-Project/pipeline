
#include <cmath>
#include <vector>
#include <iostream>
#include <cstring>
#include <fstream>
#include <images.hpp>
#include <mycomplex.hpp>
#include "dedispersion_gpu.hpp"

#define NTHREADS 1024

__device__
float compute_sweep(Complex<float>* images, unsigned long long image_size, int *delay_row, int start_time_step, int start_frequency_band,
        int freq_batch_size, int freq_batch_size_subtract, int start_band_delay, int max_steps, unsigned long long pix_idx){
    float partial_integral = 0;
    // Start the sweep from the top start frequency going downwards
    int end_frequency = start_frequency_band - freq_batch_size + 1 + freq_batch_size_subtract;
    const size_t n_pixels_in_time_int {image_size * freq_batch_size};
    if(end_frequency < 0) end_frequency = 0;
    for(int current_band_idx = start_frequency_band; current_band_idx >= end_frequency; current_band_idx--){
        int cumulative_delay = delay_row[current_band_idx + 1] - delay_row[start_frequency_band +1];
        if (current_band_idx < start_frequency_band){
            // this is needed when we start from a "sub-delay"
            cumulative_delay -= (delay_row[start_frequency_band] - delay_row[start_frequency_band +1] - start_band_delay);
        }
        // get number of (extra) timesteps required to be spent in current band
        int delay_steps = start_frequency_band == current_band_idx ? start_band_delay : (delay_row[current_band_idx] - delay_row[current_band_idx + 1]);
        // are there any steps outside the current batch window?
        int exess_steps = cumulative_delay + 1 + delay_steps - max_steps;
        if (exess_steps > 0)
            delay_steps -= exess_steps;

        // get pixels related to current freq and delay
        for(int t {0}; t <= delay_steps; t++){
            int time_step = start_time_step + cumulative_delay + t;
            Complex<float> *image_start = images + n_pixels_in_time_int * time_step + image_size * (current_band_idx - end_frequency);
            partial_integral += image_start[pix_idx].real;
        }
        if(exess_steps > 0){
            // We arrived at the end of the sweep, return 
            return partial_integral;
        }
    }
    return partial_integral;
}


__global__
void compute_partial_dedispersion_kernel(Complex<float>* images, unsigned long long side_size, int start_freq_idx, int freq_batch_size, int n_frequencies, int batch_size, int *delay_table, float *dm_starttime, 
    int n_dms, int table_size, int window_start_idx, int window_step_offset){

    unsigned int gridSize = gridDim.x * blockDim.x;

    auto cti = [window_start_idx, table_size](int offset) {
        return (window_start_idx + offset) % table_size;
    };

    long long cell_size = n_dms * table_size;
    long long width = cell_size * side_size;
    

    int end_freq = start_freq_idx - freq_batch_size + 1;
    if(end_freq < 0) end_freq = 0;

    for(int dm_idx = 0; dm_idx < n_dms; dm_idx++){
        int *delay_row = delay_table + dm_idx * n_frequencies;
        // cumulative delay before the current band to analyse
        int upper_delay = delay_row[start_freq_idx+1];
        // The following are the new sweeps
        for(int ts {0}; ts < batch_size; ts++){
            int start_time = window_step_offset + ts - upper_delay;
            if(start_time < 0) continue;
            
            int start_delay = (delay_row[start_freq_idx] - delay_row[start_freq_idx + 1]);
            unsigned long long pix_id {blockIdx.x * blockDim.x + threadIdx.x};
            unsigned long long n_pixels {side_size * side_size};
            for(; pix_id < n_pixels; pix_id += gridSize){
                int x = static_cast<int>(pix_id % side_size);
                int y = static_cast<int>(pix_id / side_size);
                float partial_integral = compute_sweep(images, n_pixels, delay_row, ts, start_freq_idx, freq_batch_size, 0, start_delay, batch_size - ts, pix_id);
                //dm_starttime[y * width + x*cell_size + dm_idx * table_size +  cti(start_time)] += partial_integral;
                atomicAdd(&dm_starttime[y * width + x*cell_size + dm_idx * table_size +  cti(start_time)], partial_integral);
            }
        }
        // continue sweeps started in previous batches
        // for each frequency band, including the top one..
        for(int start_frequency_band_idx {start_freq_idx}, a = 0; start_frequency_band_idx >= end_freq; start_frequency_band_idx--, a++){
            int delay_in_prev_bands = delay_row[ start_frequency_band_idx+1];
            // for each delay in such band that is greater than zero..
            // (if it were zero, then the sweep would have not started at a previous time step or batch)
            int delay = delay_row[ start_frequency_band_idx] - delay_row[start_frequency_band_idx + 1];
            for(int delay_offset = 1; delay_offset <= delay; delay_offset++){
                int delay_current_batch = delay_offset;
                int delay_prev_batch = delay - delay_offset;
                int start_time = window_step_offset - delay_in_prev_bands - delay_prev_batch - 1;
                if (start_time < 0) {
                    // start time is outside the window
                    // This will only happen at the start
                    continue;
                }
                // If we are at this point, then we need to sum # delay_current_batch pixels in the current
                // frequency before continuing to the following bands as usual
                unsigned long long pix_id {blockIdx.x * blockDim.x + threadIdx.x};
                unsigned long long n_pixels {side_size * side_size};
                for(; pix_id < n_pixels; pix_id += gridSize){
                    int x = static_cast<int>(pix_id % side_size);
                    int y = static_cast<int>(pix_id / side_size);
                    float partial_integral = compute_sweep(images, n_pixels, delay_row, 0, start_frequency_band_idx, freq_batch_size, a, delay_current_batch - 1, batch_size, pix_id);
                    //dm_starttime[y * width + x*cell_size + dm_idx * table_size +  cti(start_time)] += partial_integral;
                    atomicAdd(&dm_starttime[y * width + x*cell_size + dm_idx * table_size +  cti(start_time)], partial_integral);
                }
            }
        }
    }
}


void compute_partial_dedispersion_gpu(Images& images, int start_freq_idx, int freq_batch_size, int n_frequencies, int batch_size, int *delay_table, float *dm_starttime, 
        int n_dms, int table_size, int window_start_idx, int window_step_offset){
    
    struct gpuDeviceProp_t props;
    int gpu_id = -1;
    gpuGetDevice(&gpu_id);
    gpuGetDeviceProperties(&props, gpu_id);
    unsigned int n_blocks = props.multiProcessorCount * 2;
    compute_partial_dedispersion_kernel<<<n_blocks, NTHREADS>>>(reinterpret_cast<Complex<float>*>(images.data()), images.side_size, start_freq_idx, freq_batch_size,
        n_frequencies, batch_size, delay_table, dm_starttime, n_dms, table_size, window_start_idx, window_step_offset);
    gpuCheckLastError();
    gpuDeviceSynchronize();
}
