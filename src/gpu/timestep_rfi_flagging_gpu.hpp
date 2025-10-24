#ifndef __TIMESTEP_RFI_GPU__
#define __TIMESTEP_RFI_GPU__

#include <vector>
#include <images.hpp>

std::vector<float> compute_timestep_rms_gpu(Images& images);

#endif