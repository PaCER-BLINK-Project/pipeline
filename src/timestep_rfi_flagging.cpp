
#include <cmath>
#include <deque>
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


size_t flag_timestep_using_history(Images& images, std::deque<float>& history_rms, size_t history_length, unsigned int threshold){
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
    size_t count {0u};
    if(history_rms.size() < history_length){
        for(float e : rms_vector) history_rms.push_back(e);
    }else{
        std::vector<float> tmp (history_rms.size(), 0.0f);
        for(size_t i {0}; i < tmp.size(); i++) tmp[i] = history_rms[i];
        std::pair<float, float> median_rms_of_history = compute_iqr_rms(tmp);
        
        auto& flags = images.get_flags();
        if(flags.size() == 0) flags.resize(images.size(), false);

        std::vector<size_t> flagged_ts;
        float min_rms {*std::min_element(history_rms.begin(), history_rms.end())};
        for(size_t i {0}; i < rms_vector.size(); i++){
            if(rms_vector[i] > 1.3 * min_rms){
                flagged_ts.push_back(i);
                // flag all images in the timestep
                std::cout << "flag_timestep_rfi: flagging step " << i << std::endl;
                for(size_t j {0}; j < images.n_channels; j++)
                    flags[i * images.n_channels + j] = true;
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
                    flags[i * images.n_channels + j] = true;
                    
            }
        }
        for(float e : rms_vector){
            history_rms.push_back(e);
            history_rms.pop_front();
        }
        for(bool b : flags) if(b) count++;
    }
    return count;
}


/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
size_t flag_timestep_rfi(Images& images, double rms_threshold, bool use_iqr){
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
    auto& flags = images.get_flags();
    if(flags.size() == 0) flags.resize(images.size(), false);
    
    std::vector<size_t> flagged_ts;
    // now compute avg rms of all rms
    for(size_t i {0}; i < rms_vector.size(); i++){
        if(rms_vector[i] > min_rms * 3){
            flagged_ts.push_back(i);
            // flag all images in the timestep
            std::cout << "flag_timestep_rfi: flagging step " << i << std::endl;
            for(size_t j {0}; j < images.n_channels; j++)
                flags[i * images.n_channels + j] = true;
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
                flags[i * images.n_channels + j] = true;
        }
    }
    size_t count {0u};
    for(bool b : flags) if(b) count++;
    return count;
}