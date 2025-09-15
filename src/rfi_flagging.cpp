#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <limits>
#include <mycomplex.hpp>
#include <images.hpp>


template <typename T>
T compute_iqr_rms(const std::vector<T>& values){
    std::vector<T> sorted_values {values};
    std::sort(sorted_values.begin(), sorted_values.end());
    size_t count = sorted_values.size();
    size_t q75 = static_cast<size_t>(count*0.75);
    size_t q25 = static_cast<size_t>(count*0.25);
    T iqr = sorted_values[q75] - sorted_values[q25];
    return iqr / 1.35;
}


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

/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
void flag_rfi(Images& images, double rms_threshold, int radius, bool compute_iqr){
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