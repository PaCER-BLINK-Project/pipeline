
#include <cmath>
#include <iostream>
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


size_t flag_values(std::vector<float>& values, std::vector<bool>& flags, double threshold, std::string name, size_t n_intervals, size_t n_channels){
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



/**
    @brief: computes a vector of boolean values indicating which of the input images are contaminated
    by RFI. The function uses a high RMS as a signal for bad/corrupted channels.
*/
size_t flag_rfi(Images& images, double mean_threshold, double rms_threshold){
    std::vector<float> rms_vector;
    std::vector<float> mean_vector;
    
    #ifdef __GPU__
    if(gpu_support() && num_available_gpus() > 0 && images.on_gpu()){
        compute_images_mean_rms_gpu(images, mean_vector, rms_vector);
    }else{
        compute_images_mean_rms_cpu(images, mean_vector, rms_vector);
    }
    #else
    compute_images_mean_rms_cpu(images, mean_vector, rms_vector);
    #endif
    
    auto& flags = images.get_flags();
    if(flags.size() == 0) flags.resize(images.size(), false);
    
    size_t count {0u};
    count += flag_values(rms_vector, flags, rms_threshold, "rms", images.n_intervals, images.n_channels);
    std::cout << "Flagged after rms: " << count << std::endl;
    count += flag_values(mean_vector, flags, mean_threshold, "mean", images.n_intervals, images.n_channels);
    std::cout << "Flagged after mean: " << count << std::endl;
    return count;
}