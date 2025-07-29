
#include <cmath>
#include <vector>
#include <iostream>
#include <cstring>
#include <fstream>
#include <images.hpp>
#include "dedispersion.hpp"

const float K = 4.15;


float dispersive_delay_s(float DM, float f_low_ghz, float f_high_ghz){
    return K * DM * (std::pow(f_low_ghz, -2) - std::pow(f_high_ghz, -2)) / 1000.0f;
}


std::vector<int> compute_delay_table(const std::vector<float>& frequencies, const std::vector<float>& dm_list, float int_time) {
    int n_frequencies = frequencies.size();
    int n_dms = dm_list.size();
    auto delay_table = std::vector<int> (n_dms * n_frequencies, 0);
    int top_freq_idx = n_frequencies - 1;
    for(int dm_idx = 0; dm_idx < n_dms; dm_idx++){
        float dm = dm_list[dm_idx];
        delay_table[dm_idx * n_frequencies + top_freq_idx] = 0;
        for(int i = n_frequencies - 2; i >= 0; i--){
            delay_table[dm_idx * n_frequencies + i] = static_cast<int>(std::roundl((dispersive_delay_s(dm, frequencies[i], frequencies[top_freq_idx]) / int_time)));
        }
    }
    return delay_table;
}

std::vector<float> compute_normalisation_factor(std::vector<int>&delay_table, int n_dms, int n_frequencies){
    float factor = 0.0f;
    std::vector<float> norm_factors(n_dms, 0.0f);
    for(int dm_idx = 0; dm_idx < n_dms; dm_idx++){
        for(int f = 0; f < n_frequencies - 1; f++){
            norm_factors[dm_idx] += (delay_table[dm_idx * n_frequencies + f] - delay_table[dm_idx * n_frequencies + f + 1] + 1);
        }
    }
    return norm_factors;
}

float compute_sweep(Images& images, int *delay_row, int start_time_step, int start_frequency_band,
        int freq_batch_size, int start_band_delay, int max_steps, int x, int y){
    float partial_integral = 0;
    // Start the sweep from the top start frequency going downwards
    int end_frequency = start_frequency_band - freq_batch_size + 1;
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
            std::complex<float> *image_start = images.at(time_step, current_band_idx - end_frequency);
            // const int x = 59, y = 629;
            partial_integral += image_start[y * images.side_size + x].real();
        }
        
        if(exess_steps > 0){
            // We arrived at the end of the sweep, return 
            return partial_integral;
        }
    }
    return partial_integral;
}



void compute_partial_dedispersion(Images& images, int start_freq_idx, int freq_batch_size, int n_frequencies, int batch_size, int *delay_table, float *dm_starttime, 
    int n_dms, int table_size, int window_start_idx, int window_step_offset){

    auto cti = [window_start_idx, table_size](int offset) {
        return (window_start_idx + offset) % table_size;
    };

    unsigned long long cell_size = n_dms * table_size;
    unsigned long long width = cell_size * images.side_size;
    

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
            #pragma omp parallel for collapse(2)
            for(int x = 0; x < images.side_size; x++){
                for(int y = 0; y < images.side_size; y++){
                    float partial_integral = compute_sweep(images, delay_row, ts, start_freq_idx, freq_batch_size, start_delay, batch_size - ts, x, y);
                    dm_starttime[y * width + x*cell_size + dm_idx * table_size +  cti(start_time)] += partial_integral;
                }
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
                #pragma omp parallel for collapse(2)
                for(int x = 0; x < images.side_size; x++){
                    for(int y = 0; y < images.side_size; y++){
                        float partial_integral = compute_sweep(images, delay_row, 0, start_frequency_band_idx, freq_batch_size - a, delay_current_batch - 1, batch_size, x, y);
                        dm_starttime[y * width + x*cell_size + dm_idx * table_size +  cti(start_time)] += partial_integral;
                    }
                }
            }
        }
    }
}

void clear_buffer(float* dm_starttime, int side_size, int n_dms, int start_idx, int table_size, int buffer_size){
    unsigned long long  width = n_dms * table_size * side_size;
    unsigned long long  cell_size = n_dms * table_size;
    #pragma omp parallel for collapse(2)
    for(int x = 0; x < side_size; x++){
        for(int y = 0; y < side_size; y++){
            for(int dm {0}; dm < n_dms; dm++){
                for(int i {0}; i < buffer_size; i++){
                    int adj_i = (start_idx + i) % table_size;
                    dm_starttime[y * width + x*cell_size + dm * table_size + adj_i] = 0;
                }
            }
        }
    }
}

void get_elements(float *dm_starttime, int side_size, int n_dms, int dm_idx, int window_start_idx, int buffer_size, 
        int table_size, float norm_factor, int x, int y){
    
    unsigned long long  width = n_dms * table_size * side_size;
    unsigned long long  cell_size = n_dms * table_size;

    for(int i {0}; i < buffer_size; i++){
        int idx = (window_start_idx + i) % table_size;
        std::cout << "BUFFOUT (dmidx = " << dm_idx << ") " << dm_starttime[y * width + x*cell_size + dm_idx * table_size + idx] / norm_factor << std::endl;
    }
}

void dump_buffer(float *dm_starttime, int side_size, int n_dms, int window_start_idx, int buffer_size, 
        int table_size, std::string filename, int x_0, int y_0){
    unsigned long long width = n_dms * table_size * side_size;
    unsigned long long cell_size = n_dms * table_size;
    std::ofstream out_file(filename, std::ios::app | std::ios::binary);
    #define write_to_file(X) out_file.write(reinterpret_cast<char*>(&X), sizeof(X))

    if(x_0 != -1 && y_0!= -1){
        int fakedim {1};
        write_to_file(fakedim);
        write_to_file(n_dms);
        write_to_file(buffer_size);
        for(int dm {0}; dm < n_dms; dm++){
            for(int i {0}; i < buffer_size; i++){
                int idx = (window_start_idx + i) % table_size;
                float val =  dm_starttime[y_0 * width + x_0*cell_size + dm * table_size + idx];
                write_to_file(val);
            }
        }
    }else{
        // dump all pixels
        write_to_file(side_size);
        write_to_file(n_dms);
        write_to_file(buffer_size);
        for(int x = 0; x < side_size; x++){
            for(int y = 0; y < side_size; y++){
                for(int dm {0}; dm < n_dms; dm++){
                    for(int i {0}; i < buffer_size; i++){
                        int idx = (window_start_idx + i) % table_size;
                        float val =  dm_starttime[y * width + x*cell_size + dm * table_size + idx];
                        write_to_file(val);
                    }
                }
            }
        }
    }
    out_file.close();
}
