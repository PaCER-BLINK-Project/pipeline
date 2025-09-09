#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
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



float compute_image_rms(const Complex<float>* data, int side_size, int radius, bool compute_iqr){
    double sum {0.0};
    double sum2 {0.0};
    int cnt {0};
    double iqr {0.0};
    double rms_iqr {0.0};
    double median {0.0};   
    double rms {0.0}, mean {0.0};
    int center = side_size / 2;
    
    int start = center - radius;
    int end = center + radius;

    std::vector<float> iqr_values;
    int nan_count = 0, total_count = 0;

    for(int y = start; y < end; y++){
        for(int x = start; x < end; x++){
            float val = data[y * side_size + x].real;
            total_count++;
            if(std::isnan(val) || std::isinf(val)){
                nan_count++;
                continue;
            }
            if(compute_iqr) iqr_values.push_back(val);
            sum  += val;
            sum2 += val*val;
            cnt  += 1;
        }
    }
    
    mean = sum / cnt;
    rms = std::sqrt(sum2/cnt - mean*mean);

    if(compute_iqr) return compute_iqr_rms(iqr_values);
    return rms;
}

/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
std::vector<bool> flag_rfi(Images& images, double rms_threshold, int radius, bool compute_iqr){
    size_t n_images = images.size();
    std::vector<float> rms_vector (n_images);
    #pragma omp parallel for collapse(2) schedule(static)
    for(size_t i {0u}; i < images.integration_intervals(); i++){
        for(size_t j {0u}; j < images.nFrequencies; j++){
            rms_vector[i * images.nFrequencies + j] = compute_image_rms(reinterpret_cast<Complex<float>*>(images.at(i, j)), images.side_size, radius, compute_iqr);
        }
    }
    float rms = compute_iqr_rms(rms_vector);
    // print some stats
    std::cout << "flag_rfi - "
    "min rms: " << *std::min_element(rms_vector.begin(), rms_vector.end()) << ", "
    "max rms: " << *std::max_element(rms_vector.begin(), rms_vector.end()) << ", rms = " << rms <<  std::endl;
    std::vector<bool> flags (n_images, false);
    // now compute avg rms of all rms
    for(size_t i {0}; i < rms_vector.size(); i++){
        if(rms_vector[i] > rms_threshold * rms) flags[i] = true; 
    }
    return flags;
}