#ifndef __DEDISPERSION_H__
#define __DEDISPERSION_H__

#include <images.hpp>
#include <vector>

void compute_partial_dedispersion(Images& images, int start_freq_idx, int freq_batch_size, int n_frequencies, 
    int batch_size, int *delay_table, float *dm_starttime, int n_dms, int table_size, int window_start_idx, int window_step_offset);

std::vector<int> compute_delay_table(const std::vector<float>& frequencies, const std::vector<float>& dm_list, float int_time);

void clear_buffer(float* dm_starttime, int side_size, int n_dms, int start_idx, int table_size, int buffer_size);

void get_elements(float *dm_starttime, int side_size, int n_dms, int dm_idx, int window_start_idx, int buffer_size, int table_size, float norm_factor, int x, int y);

void dump_buffer(float *dm_starttime, int side_size, int n_dms, int window_start_idx, int buffer_size, 
        int table_size, const std::vector<float>& norm_factors, std::string filename, int x_0=-1, int y_0=-1);
        
std::vector<float> compute_normalisation_factor(std::vector<int>&delay_table, int n_dms, int n_frequencies);
#endif