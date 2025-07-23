#ifndef __DEDISPERSION_GPU_H__
#define __DEDISPERSION_GPU_H__

#include <images.hpp>

void compute_partial_dedispersion_gpu(Images& images, int start_freq_idx, int freq_batch_size, int n_frequencies, 
    int batch_size, int *delay_table, float *dm_starttime, int n_dms, int table_size, int window_start_idx, int window_step_offset);

#endif