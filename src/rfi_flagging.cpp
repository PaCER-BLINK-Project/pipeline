
#include <cmath>
#include <iostream>
#include <set>
#include <mycomplex.hpp>
#include <images.hpp>
#include "rfi_flagging.hpp"
#ifdef __GPU__
#include "gpu/rfi_flagging_gpu.hpp"
#endif


std::pair<float, float> compute_image_mean_rms_cpu(const Complex<float>* data, int side_size, int radius, bool compute_iqr){
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


void compute_images_mean_rms_cpu(Images& images, std::vector<float>& mean, std::vector<float>& rms){
    size_t n_images = images.size();
    mean.resize(n_images);
    rms.resize(n_images);

    #pragma omp parallel for collapse(2) schedule(static)
    for(size_t i {0u}; i < images.integration_intervals(); i++){
        for(size_t j {0u}; j < images.n_channels; j++){
            auto median_rms = compute_image_mean_rms_cpu(reinterpret_cast<Complex<float>*>(images.at(i, j)), images.side_size, 50, false);
            mean[i * images.n_channels + j] = median_rms.first;
            rms[i * images.n_channels + j] = median_rms.second;
        }
    }
}


void compute_timestep_mean_rms_cpu(std::vector<float>& images_mean, std::vector<float>& images_rms, size_t n_intervals,
        size_t n_channels, std::vector<float>& timestep_mean_of_mean, std::vector<float>& timestep_mean_of_rms){

    timestep_mean_of_mean.resize(n_intervals);
    //timestep_rms_of_mean.resize(n_intervals);
    timestep_mean_of_rms.resize(n_intervals);
    //timestep_rms_of_rms.resize(n_intervals);

    #pragma omp parallel for schedule(static)
    for(size_t i {0u}; i < n_intervals; i++){
        std::vector<float> mean_values;
        std::vector<float> rms_values;
        for(size_t j {0u}; j < n_channels; j++){
            size_t idx {i * n_channels + j};
            mean_values.push_back(images_mean[idx]);
            rms_values.push_back(images_rms[idx]);
        }
        auto p1 = compute_running_rms(mean_values, -1);
        auto p2 = compute_running_rms(rms_values, -1);

        timestep_mean_of_mean[i] = p1.first;
        //timestep_rms_of_mean[i] = p1.second;
        timestep_mean_of_rms[i] = p2.first;
        //timestep_rms_of_rms[i] = p2.second;
    }
}



void compute_channel_mean_rms_cpu(std::vector<float>& images_mean, std::vector<float>& images_rms, size_t n_intervals,
        size_t n_channels, std::vector<float>& channel_mean_of_mean, std::vector<float>& channel_mean_of_rms){

    channel_mean_of_mean.resize(n_channels);
    //timestep_rms_of_mean.resize(n_intervals);
    channel_mean_of_rms.resize(n_channels);
    //timestep_rms_of_rms.resize(n_intervals);

    #pragma omp parallel for schedule(static)
    for(size_t j {0u}; j < n_channels; j++){
        std::vector<float> mean_values;
        std::vector<float> rms_values;
        for(size_t i {0u}; i < n_intervals; i++){
            size_t idx {i * n_channels + j};
            mean_values.push_back(images_mean[idx]);
            rms_values.push_back(images_rms[idx]);
        }
        auto p1 = compute_running_rms(mean_values, -1);
        auto p2 = compute_running_rms(rms_values, -1);

        channel_mean_of_mean[j] = p1.first;
        //timestep_rms_of_mean[i] = p1.second;
        channel_mean_of_rms[j] = p2.first;
        //timestep_rms_of_rms[i] = p2.second;
    }
}



size_t flag_values(std::vector<float>& values, std::vector<bool>& flags, double threshold, std::string name){
    std::pair<float, float> median_rms_of_values = compute_iqr_rms(values);

    // print some stats
    std::cout << "flag_rfi - " << name << " - "
    "min: " << *std::min_element(values.begin(), values.end()) << ", "
    "max: " << *std::max_element(values.begin(), values.end()) << ", "
    "median: " << median_rms_of_values.first <<  ", "
    "rms: " << median_rms_of_values.second <<  std::endl;

    // now compute avg rms of all rms
    size_t count {0u};
    if(name == "rms"){
        for(size_t i {0u}; i < values.size(); i++){
            if(values[i] > median_rms_of_values.first + threshold * median_rms_of_values.second){
                flags[i] = true;
                count++;
            }
        }
    }else{
        double up {median_rms_of_values.first + threshold * median_rms_of_values.second};
        double down {median_rms_of_values.first - threshold * median_rms_of_values.second};
        for(size_t i {0u}; i < values.size(); i++){
            if(values[i] > up || values[i] < down){
                flags[i] = true;
                count++;
            }
        }
    }
    return count;
}



void flag_timestep_values(std::vector<float>& values, std::set<unsigned int>& flagged_ts, double threshold, std::string name){
    std::pair<float, float> median_rms_of_values = compute_iqr_rms(values);
    float min_value {*std::min_element(values.begin(), values.end())};
    std::cout << "flag_timestep_value - " << name << " - "
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
    std::cout << "flag_timestep_value - flagged " << flagged_ts.size() << "/" << values.size() << " of timesteps." << std::endl;
}



void flag_values_with_history(std::vector<float>& values, std::vector<float>& history, std::set<unsigned int>& flagged, double threshold){
    std::pair<float, float> median_rms_of_history = compute_iqr_rms(history);
    std::cout << "flag_timestep_rfi_with_history - "
        "min rms: " << *std::min_element(history.begin(), history.end()) << ", "
        "max rms: " << *std::max_element(history.begin(), history.end()) << ", "
        "median = " << median_rms_of_history.first << ", "
        "rms = " << median_rms_of_history.second <<  std::endl;

    for(unsigned int i {0}; i < values.size(); i++){
        if((values[i] > median_rms_of_history.first + threshold * median_rms_of_history.second) || 
                (values[i] < median_rms_of_history.first - threshold * median_rms_of_history.second)){
            flagged.insert(i);
        }
    }
    if(flagged.size() > 4){
        // flags all the steps in between the minimum and maximum flagged timestep
        // as an aggressive measure
        unsigned int min_ts {*std::min_element(flagged.begin(), flagged.end())};
        unsigned int max_ts {*std::max_element(flagged.begin(), flagged.end())};
        for(unsigned int i {min_ts}; i < max_ts; i++) flagged.insert(i);
    }
}


/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
size_t flag_rfi(Images& images, double threshold, std::deque<std::pair<float, float>>& history, int history_length){

    std::vector<float> images_rms;
    std::vector<float> images_mean;
    
    #ifdef __GPU__
    if(gpu_support() && num_available_gpus() > 0 && images.on_gpu()){
        compute_images_mean_rms_gpu(images, images_mean, images_rms);
    }else{
        compute_images_mean_rms_cpu(images, images_mean, images_rms);
    }
    #else
    compute_images_mean_rms_cpu(images, images_mean, images_rms);
    #endif
    
    auto& flags = images.get_flags();
    if(flags.size() == 0) flags.resize(images.size(), false);
    
    size_t count {0u};
    // TODO check if there are enough images a block.
    // Step 1: flag images with high variation in rms and mean values
    count += flag_values(images_rms, flags, threshold, "rms");
    std::cout << "Flagged after rms: " << count << std::endl;
    count += flag_values(images_mean, flags, threshold, "mean");
    std::cout << "Flagged after mean: " << count << std::endl;


    std::vector<float> timestep_mean_of_mean, timestep_mean_of_rms, channel_mean_of_mean, channel_mean_of_rms;

    compute_timestep_mean_rms_cpu(images_mean, images_rms, images.n_intervals, images.n_channels,
        timestep_mean_of_mean, timestep_mean_of_rms);
    
    compute_channel_mean_rms_cpu(images_mean, images_rms, images.n_intervals, images.n_channels,
        channel_mean_of_mean, channel_mean_of_rms);
    
    std::set<unsigned int> flagged_ts;
    flag_timestep_values(timestep_mean_of_mean, flagged_ts, threshold, "mean");
    flag_timestep_values(timestep_mean_of_rms, flagged_ts,  threshold, "rms");

    std::set<unsigned int> flagged_ch;
    flag_timestep_values(channel_mean_of_mean, flagged_ch, threshold, "mean");
    flag_timestep_values(channel_mean_of_rms, flagged_ch, threshold, "rms");

    /*
    if(history_length >= 0){
        if(static_cast<int>(history.size()) < history_length){
            for(size_t i {0}; i < timestep_mean_of_mean.size(); i++) {
                history.push_back({timestep_mean_of_mean[i], timestep_mean_of_rms[i]});
            }
        }else{
            std::vector<float> history_rms (history.size(), 0.0f);
            std::vector<float> history_mean (history.size(), 0.0f);
            for(size_t i {0}; i < history.size(); i++) {
                history_mean[i] = history[i].first;
                history_rms[i] = history[i].second;
            }
            flag_values_with_history(timestep_mean_of_mean, history_mean, flagged_ts, threshold);
            flag_values_with_history(timestep_mean_of_rms, history_rms, flagged_ts, threshold);
            
            for(size_t i {0}; i < timestep_mean_of_mean.size(); i++) {
                history.push_back({timestep_mean_of_mean[i], timestep_mean_of_rms[i]});
                history.pop_front();
            }
        }        
    }*/

    
    for(size_t j : flagged_ch){
        for(size_t i {0}; i < images.n_intervals; i++){
            flags[i * images.n_channels + j] = true;
            count++;
        }
    }

    for(size_t i : flagged_ts){
        for(size_t j {0}; j < images.n_channels; j++){
            flags[i * images.n_channels + j] = true;
            count++;
        }
    }
    
    return count;
}