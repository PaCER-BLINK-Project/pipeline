#ifndef __PEAKFINDING_H__
#define __PEAKFINDING_H__

#include <vector>
#include <string>

void peakfinding_simple_avg(float *dm_starttime, int side_size, std::vector<float>& dm_list, int window_start_idx, int buffer_size, 
        int table_size, int global_offset, const std::vector<float>& norm_factors, float SNR, std::string filename);

#endif