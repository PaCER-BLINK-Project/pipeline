
#include <cmath>
#include <iostream>
#include <mycomplex.hpp>
#include <images.hpp>
#include "rfi_flagging.hpp"
#ifdef __GPU__
#include "gpu/timestep_rfi_flagging_gpu.hpp"
#endif


std::vector<float> compute_timestep_rms_cpu(Images& images){
    size_t centre {images.side_size / 2};
    std::vector<float> rms_vector(images.n_intervals);
    #pragma omp parallel for schedule(static)
    for(size_t i {0u}; i < images.n_intervals; i++){
        std::vector<float> values;
        for(size_t j {0u}; j < images.n_channels; j++){
            float val {(images.at(i, j) + centre * images.side_size + centre)->real()};
            values.push_back(val);
        }
        auto p = compute_running_rms(values, -1);
        rms_vector[i] = p.second;
    }
    return rms_vector;
}


/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
void flag_timestep_rfi(Images& images, double rms_threshold, bool use_iqr){
    std::vector<float> rms_vector;
    #ifdef __GPU__
    if(gpu_support() && num_available_gpus() > 0 && images.on_gpu()){
        rms_vector = compute_timestep_rms_gpu(images);
    }else{
        rms_vector = compute_timestep_rms_cpu(images);
    }
    #else
    rms_vector = compute_timestep_rms_cpu(images);
    #endif
    
    std::pair<float, float> median_rms_of_rms = use_iqr ? compute_iqr_rms(rms_vector) : compute_running_rms(rms_vector);
    float min_rms {*std::min_element(rms_vector.begin(), rms_vector.end())};
    std::cout << "flag_timestep_rfi - "
    "min rms: " << *std::min_element(rms_vector.begin(), rms_vector.end()) << ", "
    "max rms: " << *std::max_element(rms_vector.begin(), rms_vector.end()) << ", "
    "median = " << median_rms_of_rms.first << ", "
    "rms = " << median_rms_of_rms.second <<  std::endl;
    std::vector<bool> new_flags(images.size(), false);
    auto existing_flags = images.get_flags();
    if(existing_flags.size() > 0){
        std::cout << "flag_timestep_rfi: flags are already there.." << std::endl;
        for(int i {0}; i < existing_flags.size(); i++)
            new_flags[i] = existing_flags[i];
    }
    
    std::vector<size_t> flagged_ts;
    // now compute avg rms of all rms
    for(size_t i {0}; i < rms_vector.size(); i++){
        if(rms_vector[i] > min_rms * 3){
            flagged_ts.push_back(i);
            // flag all images in the timestep
            std::cout << "flag_timestep_rfi: flagging step " << i << std::endl;
            for(size_t j {0}; j < images.n_channels; j++)
                new_flags[i * images.n_channels + j] = true;
        }
    }
    if(flagged_ts.size() > 4){
        // flags all the steps in between the minimum and maximum flagged timestep
        // as an aggressive measure
        size_t min_ts {*std::min_element(flagged_ts.begin(), flagged_ts.end())};
        size_t max_ts {*std::max_element(flagged_ts.begin(), flagged_ts.end())};
        for(size_t i {min_ts}; i < max_ts; i++){
            // flag all images in the timestep
            std::cout << "flag_timestep_rfi: extra flagging step " << i << std::endl;
            for(size_t j {0}; j < images.n_channels; j++)
                new_flags[i * images.n_channels + j] = true;
        }
    }
    images.set_flags(new_flags);
}