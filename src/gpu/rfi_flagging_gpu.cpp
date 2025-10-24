#include <vector>
#include <cmath>
#include <iostream>
#include <mycomplex.hpp>
#include <images.hpp>
#include <math.h>
#include <memory_buffer.hpp>

#define CHECK_RADIUS 50
#define SELECTION_SIDE (CHECK_RADIUS * 2 + 1)
#define RMS_THRESHOLD 4

#define NTHREADS 1024
#ifdef __NVCC__
#define WARPSIZE 32
#else
#define WARPSIZE 64
#endif
#define NWARPS (NTHREADS / WARPSIZE)




__global__ void compute_images_rms_kernel(const Complex<float>*data, unsigned int side_size, unsigned int n_images, float *out_rms){

    const unsigned int grid_size {blockDim.x * gridDim.x};
    const unsigned int lane_id {threadIdx.x % warpSize};
    const unsigned int warp_id {threadIdx.x / warpSize};
    const unsigned int image_size {side_size * side_size};
    const unsigned int center {side_size / 2};
    const unsigned int start_offset {center - CHECK_RADIUS};
    /* 
    Thread blocks are assigned disjoint sets of images to analyse.
    For each image assigned to a block central pixels are copied in 
    shared memory so that RMS can be computed over them.
    */
    __shared__ float centre_pixels[SELECTION_SIDE * SELECTION_SIDE];
    // support array to compute block-wide reduction
    __shared__ float support[3 * NWARPS];
    
    for(size_t image_id {blockIdx.x}; image_id < n_images; image_id += gridDim.x){
        const Complex<float>* image {data + image_id * image_size};
        // ==============================================================================
        // STEP 1: copy data in shared memory, but also perform a a first step reduction.
        // ==============================================================================
        float sum {0.0f}, sum2 {0.0f};
        for(unsigned int i {threadIdx.x}; i < SELECTION_SIDE * SELECTION_SIDE; i += blockDim.x){
            unsigned int y {start_offset + (i / SELECTION_SIDE)};
            unsigned int x {start_offset + (i % SELECTION_SIDE)};
            float val {image[y * side_size + x].real};
            sum += val;
            sum2 += val*val;
            centre_pixels[i] = val;
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
                support[0] = sum;
                support[1] = sum2;
            }
        }
        __syncthreads();
        float mean {support[0] / (SELECTION_SIDE * SELECTION_SIDE)};
        float stdev {sqrt(support[1] / (SELECTION_SIDE * SELECTION_SIDE) - mean*mean)};
        // ==============================================================================
        // STEP 2 Perform the second reduction, this time excluding outliers to compute
        // a significant stdev.
        // ==============================================================================
        sum = 0.0f; sum2 = 0.0f;
        unsigned int count {0};
        for(unsigned int i {threadIdx.x}; i < SELECTION_SIDE * SELECTION_SIDE; i += blockDim.x){
            unsigned int y {start_offset + (i / SELECTION_SIDE)};
            unsigned int x {start_offset + (i % SELECTION_SIDE)};
            float val = centre_pixels[i];
            if((std::abs(val - mean) / stdev) < RMS_THRESHOLD){
                sum += val;
                sum2 += val*val;
                count += 1;
            }
        }
        // Step 2: compute first pass RMS estimation with warp-wide reduction
        for(unsigned int i {warpSize / 2}; i >= 1; i /= 2){
            float up_sum = __gpu_shfl_down(sum, i);
            float up_sum2 = __gpu_shfl_down(sum2, i);
            unsigned int up_count = __gpu_shfl_down(count, i);
            if(lane_id < i){
                sum += up_sum;
                sum2 += up_sum2;
                count += up_count;
            }
            #ifdef __NVCC__
            __syncwarp();
            #endif
        }
        if(lane_id == 0 && warp_id > 0) {
            support[warp_id] = sum;
            support[NWARPS + warp_id] = sum2;
            support[2*NWARPS + warp_id] = count;
        }
        __syncthreads();
        if(warp_id == 0){
            if(lane_id > 0 && lane_id < NWARPS) {
                sum = support[lane_id];
                sum2 = support[NWARPS + lane_id];
                count = support[2*NWARPS + lane_id];
            }
            #ifdef __NVCC__
            __syncwarp();
            #endif
            for(unsigned int i {NWARPS / 2}; i >= 1; i /= 2){
                float up_sum = __gpu_shfl_down(sum, i, NWARPS);
                float up_sum2 = __gpu_shfl_down(sum2, i, NWARPS);
                unsigned int up_count = __gpu_shfl_down(count, i, NWARPS);
                if(lane_id < i){
                    sum += up_sum;
                    sum2 += up_sum2;
                    count += up_count;
                }
                #ifdef __NVCC__
                __syncwarp();
                #endif
            }
            if(lane_id == 0){
                float final_mean {sum / count};
                float final_stdev {sqrt(sum2 / count - final_mean*final_mean)};
                out_rms[image_id] = final_stdev;
            }
        }
        __syncthreads();
    }
}



__global__ void clear_flagged_images_kernel(Complex<float>*data, unsigned int side_size, unsigned int n_images, bool *flagged){
    const unsigned int grid_size {blockDim.x * gridDim.x};
    const unsigned int image_size {side_size * side_size};

    for(size_t image_id {blockIdx.x}; image_id < n_images; image_id += gridDim.x){
        if(!flagged[image_id]) continue;

        Complex<float>* image {data + image_id * image_size};
        for(unsigned int i {threadIdx.x}; i < image_size; i += blockDim.x){
            image[i].real = 0.0f;
            // image value is not used, so avoid the assignment operation.
        }
    }
}



/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
std::vector<float> compute_images_rms_gpu(Images& images){
    images.to_gpu();
    size_t n_images = images.size();
    MemoryBuffer<float> rms_vector_mb {n_images, MemoryType::DEVICE};
    std::vector<float> rms_vector(n_images);
    struct gpuDeviceProp_t props;
    int gpu_id = -1;
    gpuGetDevice(&gpu_id);
    gpuGetDeviceProperties(&props, gpu_id);
    unsigned int n_blocks = props.multiProcessorCount * 2;
    compute_images_rms_kernel<<<n_blocks, NTHREADS>>>(reinterpret_cast<Complex<float>*>(images.data()),
        static_cast<unsigned int>(images.side_size), static_cast<unsigned int>(n_images), rms_vector_mb.data());
    gpuCheckLastError();
    gpuMemcpy(rms_vector.data(), rms_vector_mb.data(), sizeof(float) * n_images, gpuMemcpyDeviceToHost);
    gpuDeviceSynchronize();
    return rms_vector;
}



void clear_flagged_images_gpu(Images& images){
    images.to_gpu();
    size_t n_images = images.size();
    auto& flags = images.get_flags();
    MemoryBuffer<bool> flagged_images {n_images};
    for(int i {0}; i < n_images; ++i){
        flagged_images[i] = flags[i];
        flags[i] = false;
    }
    flagged_images.to_gpu();
    struct gpuDeviceProp_t props;
    int gpu_id = -1;
    gpuGetDevice(&gpu_id);
    gpuGetDeviceProperties(&props, gpu_id);
    unsigned int n_blocks = props.multiProcessorCount * 2;
    clear_flagged_images_kernel<<<n_blocks, NTHREADS>>>(reinterpret_cast<Complex<float>*>(images.data()),
        static_cast<unsigned int>(images.side_size), static_cast<unsigned int>(n_images), flagged_images.data());
    gpuCheckLastError();
    gpuDeviceSynchronize();

}