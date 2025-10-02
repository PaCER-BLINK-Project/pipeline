#ifndef __RFI_FLAGGING_H__
#define __RFI_FLAGGING_H__

#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <images.hpp>


void flag_rfi(Images& images, double rms_threshold, bool on_gpu);


void flag_rfi_cpu(Images& images, double rms_threshold, int radius, bool compute_iqr);




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

#endif