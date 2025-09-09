#ifndef __RFI_FLAGGING_H__
#define __RFI_FLAGGING_H__

#include <vector>
#include <images.hpp>

std::vector<bool> flag_rfi(Images& images, double rms_threshold, int radius, bool compute_iqr);

#endif