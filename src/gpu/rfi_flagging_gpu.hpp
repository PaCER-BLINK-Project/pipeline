#ifndef __RFI_FLAGGING_GPU_H__
#define __RFI_FLAGGING_GPU_H__

#include <images.hpp>

std::vector<float> compute_images_rms_gpu(Images& images);

/*
    Set all flagged images to 0s.
*/
void clear_flagged_images_gpu(Images& images, MemoryBuffer<float>& good_pixels);

void update_good_pixels_gpu(Images& images, MemoryBuffer<float>& good_pixels);

#endif