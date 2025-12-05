
#include <cmath>
#include <set>
#include <deque>
#include <iostream>
#include <mycomplex.hpp>
#include <images.hpp>

#include "rfi_flagging.hpp"
#ifdef __GPU__
#include "gpu/timestep_rfi_flagging_gpu.hpp"
#endif


void compute_timestep_rms_cpu(Images& images, std::vector<float>& mean, std::vector<float>& rms){
    size_t centre {images.side_size / 2};
    rms.resize(images.n_intervals);
    mean.resize(images.n_intervals);
    #pragma omp parallel for schedule(static)
    for(size_t i {0u}; i < images.n_intervals; i++){
        std::vector<float> values;
        for(size_t j {0u}; j < images.n_channels; j++){
            float val {(images.at(i, j) + centre * images.side_size + centre)->real()};
            values.push_back(val);
        }
        auto p = compute_running_rms(values, -1);
        mean[i] = p.first;
        rms[i] = p.second;
    }
}



void flag_values(std::vector<float>& values, std::set<unsigned int>& flagged_ts, std::string name){
    std::pair<float, float> median_rms_of_values = compute_iqr_rms(values);
    float min_value {*std::min_element(values.begin(), values.end())};
    std::cout << "flag_timestep_rfi - " << name << " - "
    "min: " << *std::min_element(values.begin(), values.end()) << ", "
    "max: " << *std::max_element(values.begin(), values.end()) << ", "
    "median = " << median_rms_of_values.first << ", "
    "rms = " << median_rms_of_values.second <<  std::endl;
    
    if(name == "rms"){
        for(unsigned int i {0}; i < values.size(); i++){
            if(values[i] > min_value * 3)
                flagged_ts.insert(i);
        }
    }else{
        const size_t threshold {5u};
        for(unsigned int i {0}; i < values.size(); i++){
            if((values[i] > median_rms_of_values.first + threshold * median_rms_of_values.second) ||
                    (values[i] < median_rms_of_values.first - threshold * median_rms_of_values.second))
                flagged_ts.insert(i);
        }
    }
    
    if(flagged_ts.size() > 4){
        // flags all the steps in between the minimum and maximum flagged timestep
        // as an aggressive measure
        unsigned int min_ts {*std::min_element(flagged_ts.begin(), flagged_ts.end())};
        unsigned int max_ts {*std::max_element(flagged_ts.begin(), flagged_ts.end())};
        for(unsigned int i {min_ts}; i < max_ts; i++) 
            flagged_ts.insert(i);
    }
}


void flag_values_with_history(std::vector<float>& values, std::vector<float>& history, std::set<unsigned int>& flagged_ts, unsigned int threshold){
    std::pair<float, float> median_rms_of_history = compute_iqr_rms(history);
    std::cout << "flag_timestep_rfi_with_history - "
        "min rms: " << *std::min_element(history.begin(), history.end()) << ", "
        "max rms: " << *std::max_element(history.begin(), history.end()) << ", "
        "median = " << median_rms_of_history.first << ", "
        "rms = " << median_rms_of_history.second <<  std::endl;

    for(unsigned int i {0}; i < values.size(); i++){
        if((values[i] > median_rms_of_history.first + threshold * median_rms_of_history.second) || 
                (values[i] < median_rms_of_history.first - threshold * median_rms_of_history.second)){
            flagged_ts.insert(i);
        }
    }
    if(flagged_ts.size() > 4){
        // flags all the steps in between the minimum and maximum flagged timestep
        // as an aggressive measure
        unsigned int min_ts {*std::min_element(flagged_ts.begin(), flagged_ts.end())};
        unsigned int max_ts {*std::max_element(flagged_ts.begin(), flagged_ts.end())};
        for(unsigned int i {min_ts}; i < max_ts; i++) flagged_ts.insert(i);
    }
}

/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
size_t flag_timestep_rfi(Images& images, double rms_threshold, std::deque<std::pair<float, float>>& history, int history_length, unsigned int threshold){
    std::vector<float> rms_timestep;
    std::vector<float> mean_timestep;
    std::vector<float> rms_channel;
    std::vector<float> mean_channel;

    #ifdef __GPU__
    if(gpu_support() && num_available_gpus() > 0 && images.on_gpu()){
        compute_timestep_channel_rms_gpu(images, mean_timestep, rms_timestep, mean_channel, rms_channel);
    }else{
        compute_timestep_channel_rms_gpu(images, mean_timestep, rms_timestep, mean_channel, rms_channel);
    }
    #else
    compute_timestep_channel_rms_gpu(images, mean_timestep, rms_timestep, mean_channel, rms_channel);
    #endif
    auto& flags = images.get_flags();
    if(flags.size() == 0) flags.resize(images.size(), false);

    std::set<unsigned int> flagged_ts;
    flag_values(rms_timestep, flagged_ts, "rms");
    flag_values(mean_timestep, flagged_ts, "mean");
    
    // let's take care of narrowband values - won't be as aggressive as "flag_balues"
    std::set<unsigned int> flagged_ch;
    // std::pair<float, float> median_rms_of_mean_channel = compute_iqr_rms(mean_channel);
    // for(unsigned int i {0}; i < mean_channel.size(); i++){
    //     if((mean_channel[i] > median_rms_of_mean_channel.first + threshold * median_rms_of_mean_channel.second) || 
    //             (mean_channel[i] < median_rms_of_mean_channel.first - threshold * median_rms_of_mean_channel.second)){
    //         flagged_ch.insert(i);
    //     }
    // }
    //flag_values(rms_channel, flagged_ch, "rms channel");
    // flag_values(mean_channel, flagged_ch, "mean channel");

    if(history_length >= 0){
        if(static_cast<int>(history.size()) < history_length){
            for(size_t i {0}; i < rms_timestep.size(); i++) {
                history.push_back({mean_timestep[i], rms_timestep[i]});
            }
        }else{
            std::vector<float> history_rms (history.size(), 0.0f);
            std::vector<float> history_mean (history.size(), 0.0f);
            for(size_t i {0}; i < history.size(); i++) {
                history_mean[i] = history[i].first;
                history_rms[i] = history[i].second;
            }
            flag_values_with_history(rms_timestep, history_rms, flagged_ts, threshold);
            flag_values_with_history(mean_timestep, history_mean, flagged_ts, threshold);
            
            for(size_t i {0}; i < rms_timestep.size(); i++) {
                history.push_back({mean_timestep[i], rms_timestep[i]});
                history.pop_front();
            }
        }        
    }

    size_t count {0u};
    for(size_t i : flagged_ts){
        for(size_t j {0}; j < images.n_channels; j++){
            flags[i * images.n_channels + j] = true;
            count++;
        }
    }
    for(size_t j : flagged_ch){
        for(size_t i {0}; i < images.n_intervals; i++){
            flags[i * images.n_channels + j] = true;
            count++;
        }
    }
    return count;
}