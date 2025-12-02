#ifndef __TIMESTEP_RFI__
#define __TIMESTEP_RFI__

size_t flag_timestep_rfi(Images& images, double rms_threshold, bool use_iqr);
size_t flag_timestep_using_history(Images& images, std::deque<float>& history_rms, size_t history_length, unsigned int threshold = 5);

#endif