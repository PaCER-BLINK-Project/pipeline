
#include <cmath>
#include <iostream>
#include <mycomplex.hpp>
#include <images.hpp>
#include "rfi_flagging.hpp"



float compute_image_rms(const Complex<float>* data, int side_size, int radius, bool compute_iqr){
    int center = side_size / 2;
    int start = center - radius;
    int end = center + radius;
    // TODO: check boundary conditions
    std::vector<float> values;

    for(int y = start; y <= end; y++){
        for(int x = start; x <= end; x++){
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
void flag_rfi_cpu(Images& images, double rms_threshold, int radius, bool compute_iqr){
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
        std::cout << "RMS VALUE: " << rms_vector[i] << std::endl;
        if(rms_vector[i] > rms_threshold * rms) flags[i] = true;
    }
    images.set_flags(flags);
}