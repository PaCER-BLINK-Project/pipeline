#ifndef __RFI_FLAGGING_GPU_H__
#define __RFI_FLAGGING_GPU_H__

#include <images.hpp>

void flag_rfi_gpu(Images& images, double rms_threshold);

/*
    Set all flagged images to 0s.
*/
void clear_flagged_images_gpu(Images& images);

#endif