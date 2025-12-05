#ifndef __TIMESTEP_RFI__
#define __TIMESTEP_RFI__
#include <deque>
#include <images.hpp>

size_t flag_timestep_rfi(Images& images, double rms_threshold, std::deque<std::pair<float, float>>& history_rms, int history_length, unsigned int threshold = 7);

#endif