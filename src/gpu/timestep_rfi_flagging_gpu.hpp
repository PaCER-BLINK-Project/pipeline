#ifndef __TIMESTEP_RFI_GPU__
#define __TIMESTEP_RFI_GPU__

#include <vector>
#include <images.hpp>

void compute_timestep_channel_rms_gpu(Images& images, std::vector<float>& mean_timestep, std::vector<float>& rms_timestep,
        std::vector<float>& mean_channel, std::vector<float>& rms_channel);

#endif