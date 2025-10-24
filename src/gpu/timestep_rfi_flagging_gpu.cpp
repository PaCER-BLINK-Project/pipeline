#include <vector>
#include <cmath>
#include <iostream>
#include <mycomplex.hpp>
#include <images.hpp>
#include <math.h>
#include <memory_buffer.hpp>



#define NTHREADS 1024
#ifdef __NVCC__
#define WARPSIZE 32
#else
#define WARPSIZE 64
#endif
#define NWARPS (NTHREADS / WARPSIZE)


// test with blink-submit.py --centre  82.203949,21.861910 --obsid 1259685792   --dir-postfix "small_oversampling" --partition gpu-dev  --imgsize 256 --dyspec 1,1 --offset 186 --duration 30 --time 01:00:00

__global__ void compute_timestep_rms_kernel(const Complex<float>*data, unsigned int side_size, unsigned int n_channels, unsigned int n_intervals, float *out_rms){

    const unsigned int lane_id {threadIdx.x % warpSize};
    const unsigned int warp_id {threadIdx.x / warpSize};
    const unsigned int image_size {side_size * side_size};
    const unsigned int n_images {n_channels * n_intervals};
    const unsigned int center {side_size / 2};

    // support array to compute block-wide reduction
    __shared__ float support[2 * NWARPS];
    
    // each block handles a time step and all channels in it
    for(unsigned int timestep {blockIdx.x}; timestep < n_intervals; timestep += gridDim.x){
        float sum {0.0f}, sum2 {0.0f};
        for(size_t channel {threadIdx.x}; channel < n_channels; channel += blockDim.x){
            size_t image_id {timestep * n_channels + channel};
            const Complex<float>* center_image {data + image_id * image_size + center * side_size + center};
            float val {center_image->real};
            sum += val;
            sum2 += val*val;
        }
        for(unsigned int i {warpSize / 2}; i >= 1; i /= 2){
            float up_sum = __gpu_shfl_down(sum, i);
            float up_sum2 = __gpu_shfl_down(sum2, i);
            if(lane_id < i){
                sum += up_sum;
                sum2 += up_sum2;
            }
            #ifdef __NVCC__
            __syncwarp();
            #endif
        }
        if(lane_id == 0 && warp_id > 0) {
            support[warp_id] = sum;
            support[NWARPS + warp_id] = sum2;
        }
        __syncthreads();
        if(warp_id == 0){
            if(lane_id > 0 && lane_id < NWARPS) {
                sum = support[lane_id];
                sum2 = support[NWARPS + lane_id];
            }
            #ifdef __NVCC__
            __syncwarp();
            #endif
            for(unsigned int i {NWARPS / 2}; i >= 1; i /= 2){
                float up_sum = __gpu_shfl_down(sum, i, NWARPS);
                float up_sum2 = __gpu_shfl_down(sum2, i, NWARPS);
                if(lane_id < i){
                    sum += up_sum;
                    sum2 += up_sum2;
                }
                #ifdef __NVCC__
                __syncwarp();
                #endif
            }
            if(lane_id == 0){
                float mean {sum / n_channels};
                float stdev {sqrt(sum2 / n_channels- mean*mean)};
                out_rms[timestep] = stdev;
            }
        }
        __syncthreads();
    }
}



/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
std::vector<float> compute_timestep_rms_gpu(Images& images){
    std::cout << "compute_timestep_rms_gpu.." << std::endl;
    images.to_gpu();
    size_t n_images = images.size();
    MemoryBuffer<float> rms_vector_mb {images.n_intervals, MemoryType::DEVICE};
    std::vector<float> rms_vector(images.n_intervals);
    struct gpuDeviceProp_t props;
    int gpu_id = -1;
    gpuGetDevice(&gpu_id);
    gpuGetDeviceProperties(&props, gpu_id);
    unsigned int n_blocks = props.multiProcessorCount * 2;
    compute_timestep_rms_kernel<<<n_blocks, NTHREADS>>>(reinterpret_cast<Complex<float>*>(images.data()),
        static_cast<unsigned int>(images.side_size), static_cast<unsigned int>(images.n_channels),
        static_cast<unsigned int>(images.n_intervals), rms_vector_mb.data());
    gpuCheckLastError();
    gpuMemcpy(rms_vector.data(), rms_vector_mb.data(), sizeof(float) * images.n_intervals, gpuMemcpyDeviceToHost);
    gpuDeviceSynchronize();
    return rms_vector;
}