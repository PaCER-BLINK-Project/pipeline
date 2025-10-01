#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <limits>
#include <mycomplex.hpp>
#include <images.hpp>
#include <math.h>

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


template <typename T>
T compute_running_rms(const std::vector<T>& values, int threshold = 4){
    size_t count_all = values.size();
    T sum_all = std::accumulate(values.begin(), values.end(), 0);
    T mean_all {sum_all / count_all};
    
    // Now calculate the variance
    auto variance_func = [&mean_all, &count_all](T accumulator, const T& val) {
        return accumulator + ((val - mean_all)*(val - mean_all) / (count_all - 1));
    };
    T variance = std::accumulate(values.begin(), values.end(), 0.0, variance_func);
    T stdev = std::sqrt(variance);

    size_t count {0};
    T sum {0}, sum2 {0}, mean {0}, rms {0};
    for(const T& val : values){
        if(!std::isnan(val) && !std::isinf(val) && (std::abs(val - mean_all) / stdev) < threshold ){
            sum += val;
            sum2 += val*val;
            count += 1;   
        }
    }
    if(count == 0) return std::numeric_limits<T>::infinity();
    mean = sum / count;
    rms = std::sqrt(sum2 / count - mean*mean);
    return rms;
}


float compute_image_rms(const Complex<float>* data, int side_size, int radius, bool compute_iqr){
    int center = side_size / 2;
    int start = center - radius;
    int end = center + radius;

    std::vector<float> values;

    for(int y = start; y < end; y++){
        for(int x = start; x < end; x++){
            float val = data[y * side_size + x].real;
            values.push_back(val);
        }
    }
    if(compute_iqr) return compute_iqr_rms(values);
    else return compute_running_rms(values);
}


__global__ void compute_images_rms(const Complex<float>*data, size_t side_size, size_t n_images, float *out_rms){

    const unsigned int grid_size {blockDim.x * gridDim.x};
    const unsigned int lane_id {threadIdx.x % warpSize};
    const unsigned int warp_id {threadIdx.x / warpSize};
    const size_t image_size {side_size * side_size};
    const int center {side_size / 2};
    const int start_offset {center - radius};
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
        for(int i {threadIdx.x}; i < SELECTION_SIDE * SELECTION_SIDE; i += blockDim.x){
            int y {start_offset + (i / SELECTION_SIDE)};
            int x {start_offset + (i % SELECTION_SIDE)};
            float val {image[y * side_size + x].real};
            sum += val;
            sum2 += val*val;
            centre_pixels[i] = val;
        }
        for(int i {warpSize / 2}; i >= 1; i /= 2){
            float up_sum = __gpu_shfl_down(sum, i);
            float up_sum2 = __gpu_shfl_down(sum2, i);
            if(laneId < i){
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
            for(int i {NWARPS / 2}; i >= 1; i /= 2){
                float up_sum = __gpu_shfl_down(sum, i, NWARPS);
                float up_sum2 = __gpu_shfl_down(sum2, i, NWARPS);
                if(laneId < i){
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
        float stdev {sqrt(support[1] / count - mean*mean)};
        // ==============================================================================
        // STEP 2 Perform the second reduction, this time excluding outliers to compute
        // a significant stdev.
        // ==============================================================================
        sum = 0.0f; sum2 = 0.0f;
        int count {0};
        for(int i {threadIdx.x}; i < SELECTION_SIDE * SELECTION_SIDE; i += blockDim.x){
            int y {start_offset + (i / SELECTION_SIDE)};
            int x {start_offset + (i % SELECTION_SIDE)};
            float val = centre_pixels[i];
            if((std::abs(val - mean) / stdev) < RMS_THRESHOLD){
                sum += val;
                sum2 += val*val;
                count += 1;
            }
        }
        // Step 2: compute first pass RMS estimation with warp-wide reduction
        for(int i {warpSize / 2}; i >= 1; i /= 2){
            float up_sum = __gpu_shfl_down(sum, i);
            float up_sum2 = __gpu_shfl_down(sum2, i);
            int up_count = __gpu_shfl_down(count, i);
            if(laneId < i){
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
            for(int i {NWARPS / 2}; i >= 1; i /= 2){
                float up_sum = __gpu_shfl_down(sum, i, NWARPS);
                float up_sum2 = __gpu_shfl_down(sum2, i, NWARPS);
                int up_count = __gpu_shfl_down(count, i, NWARPS);
                if(laneId < i){
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

/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
void flag_rfi_gpu(Images& images, double rms_threshold, int radius){
    size_t n_images = images.size();
    std::vector<float> rms_vector (n_images);
    #pragma omp parallel for collapse(2) schedule(static)
    for(size_t i {0u}; i < images.integration_intervals(); i++){
        for(size_t j {0u}; j < images.n_channels; j++){
            rms_vector[i * images.n_channels + j] = compute_image_rms(reinterpret_cast<Complex<float>*>(images.at(i, j)), images.side_size, radius, compute_iqr);
        }
    }
    float rms = compute_iqr ? compute_iqr_rms(rms_vector) : compute_running_rms(rms_vector);
    // print some stats
    std::cout << "flag_rfi - "
    "min rms: " << *std::min_element(rms_vector.begin(), rms_vector.end()) << ", "
    "max rms: " << *std::max_element(rms_vector.begin(), rms_vector.end()) << ", rms = " << rms <<  std::endl;
    std::vector<bool> flags (n_images, false);
    // now compute avg rms of all rms
    for(size_t i {0}; i < rms_vector.size(); i++){
        if(rms_vector[i] > rms_threshold * rms) flags[i] = true;
    }
    images.set_flags(flags);
}