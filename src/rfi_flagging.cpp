
#include <cmath>
#include <iostream>
#include <mycomplex.hpp>
#include <images.hpp>
#include "rfi_flagging.hpp"
#ifdef __GPU__
#include "gpu/rfi_flagging_gpu.hpp"
#endif


std::pair<float, float> compute_image_rms_cpu(const Complex<float>* data, int side_size, int radius, bool compute_iqr){
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


std::vector<float> compute_images_rms_cpu(Images& images){
    size_t n_images = images.size();
    std::vector<float> rms_vector(n_images);

    #pragma omp parallel for collapse(2) schedule(static)
    for(size_t i {0u}; i < images.integration_intervals(); i++){
        for(size_t j {0u}; j < images.n_channels; j++){
            auto median_rms = compute_image_rms_cpu(reinterpret_cast<Complex<float>*>(images.at(i, j)), images.side_size, 50, false);
            rms_vector[i * images.n_channels + j] = median_rms.second;
        }
    }
    return rms_vector;
}

/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
size_t flag_rfi(Images& images, double rms_threshold, bool use_iqr){
    std::vector<float> rms_vector;
    #ifdef __GPU__
    if(gpu_support() && num_available_gpus() > 0 && images.on_gpu()){
        rms_vector = compute_images_rms_gpu(images);
    }else{
        rms_vector = compute_images_rms_cpu(images);
    }
    #else
    rms_vector = compute_images_rms_cpu(images);
    #endif
    
    std::pair<float, float> median_rms_of_rms = use_iqr ? compute_iqr_rms(rms_vector) : compute_running_rms(rms_vector);

    // print some stats
    std::cout << "flag_rfi - "
    "min rms: " << *std::min_element(rms_vector.begin(), rms_vector.end()) << ", "
    "max rms: " << *std::max_element(rms_vector.begin(), rms_vector.end()) << ", rms = " << median_rms_of_rms.second <<  std::endl;
    auto& flags = images.get_flags();
    if(flags.size() == 0) flags.resize(images.size(), false);
    // now compute avg rms of all rms
    size_t count {0u};
    for(size_t i {0u}; i < rms_vector.size(); i++){
        if(rms_vector[i] > median_rms_of_rms.first + rms_threshold * median_rms_of_rms.second){
            flags[i] = true;
            count++;
        }
    }
    return count;
}