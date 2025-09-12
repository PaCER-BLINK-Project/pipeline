#include <stdexcept>
#include <memory_buffer.hpp>
#include <images.hpp>
#include <complex>
#include "dynamic_spectrum.hpp"


void DynamicSpectrum::add_images(const Images& images){
    if(current_offset + images.integration_intervals() > n_timesteps)
        throw std::runtime_error {"DynamicSpectrum::add_images: tried to add more points than accounted for in dynamic spectrum."};

    #pragma omp parallel for collapse(2) schedule (static)
    for(size_t interval {0}; interval < images.integration_intervals(); interval++){
        for(size_t fine_channel {0}; fine_channel < images.nFrequencies; fine_channel++){
            if(images.is_flagged(interval, fine_channel)) continue;
            const std::complex<float>* image = images.at(interval, fine_channel);
            float val = image[y * images.side_size + x].real();
            // TODO make sure lower fine channels are shown at the bottom
            this->data()[(images.obsInfo.coarse_channel_index * images.nFrequencies + fine_channel) * n_timesteps + current_offset + interval] = val;
        }
    }
}

void DynamicSpectrum::to_fits_file(std::string filename){
    FITS fitsImage;
    FITS::HDU hdu;
    hdu.set_image(data(), n_timesteps, n_frequencies);    
    fitsImage.add_HDU(hdu);
    fitsImage.to_file(filename);
}
