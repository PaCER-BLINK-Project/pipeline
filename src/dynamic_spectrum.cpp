#include <stdexcept>
#include <memory_buffer.hpp>
#include <images.hpp>
#include <complex>
#include <mycomplex.hpp>
#include "dynamic_spectrum.hpp"

#define NTHREADS 1024

namespace {

    __global__
    void add_images_kernel(const Complex<float>* images, bool* flags, size_t image_side, size_t n_images, size_t n_channels,
            int coarse_idx, size_t current_offset, size_t n_timesteps, int x, int y, float *ds){

        const unsigned int grid_size {blockDim.x * gridDim.x};
        size_t image_size {image_side * image_side};

        for(size_t image_id {blockDim.x * blockIdx.x + threadIdx.x}; image_id < n_images; image_id += grid_size){
            if(flags[image_id]) continue;
            size_t fine_channel {image_id % n_channels};
            size_t interval {image_id / n_channels};
            const Complex<float>* image_px {images + image_id * image_size + y * image_side + x};
            float val = (*image_px).real;
            ds[(coarse_idx * n_channels + fine_channel) * n_timesteps + current_offset + interval] = val;
        }
    }


    void add_images_cpu(Images& images, size_t current_offset, size_t n_timesteps, int x, int y, float* ds){
        #pragma omp parallel for collapse(2) schedule (static)
        for(size_t interval {0}; interval < images.integration_intervals(); interval++){
            for(size_t fine_channel {0}; fine_channel < images.n_channels; fine_channel++){
                if(images.is_flagged(interval, fine_channel)) continue;
                const std::complex<float>* image = images.at(interval, fine_channel);
                float val = image[y * images.side_size + x].real();
                ds[(images.obsInfo.coarse_channel_index * images.n_channels + fine_channel) * n_timesteps + current_offset + interval] = val;
            }
        }
    }


    void add_images_gpu(Images& images, size_t current_offset, size_t n_timesteps, int x, int y, float* ds){
        size_t n_images {images.size()};
        auto flags = images.get_flags();
        MemoryBuffer<bool> flagged_images {n_images};
        if(flags.size()){
            for(size_t i {0}; i < n_images; ++i) flagged_images[i] = flags[i];
        }else{
            for(size_t i {0}; i < n_images; ++i) flagged_images[i] = false;
        }
        flagged_images.to_gpu();
        struct gpuDeviceProp_t props;
        int gpu_id = -1;
        gpuGetDevice(&gpu_id);
        gpuGetDeviceProperties(&props, gpu_id);
        unsigned int n_blocks = props.multiProcessorCount * 2;
        add_images_kernel<<<n_blocks, NTHREADS>>>(reinterpret_cast<Complex<float>*>(images.data()), flagged_images.data(),
            images.side_size, n_images, images.n_channels, images.obsInfo.coarse_channel_index, current_offset, n_timesteps, x, y, ds);
        gpuCheckLastError();
        gpuDeviceSynchronize();
    }
}


void DynamicSpectrum::add_images(Images& images){
    if(current_offset + images.integration_intervals() > n_timesteps)
        throw std::runtime_error {"DynamicSpectrum::add_images: tried to add more points than accounted for in dynamic spectrum."};

    #ifdef __GPU__
    if(images.on_gpu() && gpu_support() && num_available_gpus() > 0){
        add_images_gpu(images, current_offset, n_timesteps, x, y, this->data());
    }else{
        add_images_cpu(images, current_offset, n_timesteps, x, y, this->data());
    }
    #else
    add_images_cpu(images, current_offset, n_timesteps, x, y, this->data());
    #endif
}


void DynamicSpectrum::to_fits_file(std::string filename){
    FITS fitsImage {filename, FITS::Mode::WRITE};
    FITS::HDU hdu;
    hdu.set_image(data(), n_timesteps, n_channels);
    // hdu.set_dynspec_wcs_keywords(freq_start, delta_freq, 0.00, delta_time);
    
    // add keywords specific to dynamic spectrum:
    std::string empty_string;   
    hdu.add_keyword(std::string("CTYPE2"),std::string("FREQUENCY"),empty_string);
    hdu.add_keyword(std::string("CUNIT2"),std::string("MHz"),empty_string);
    hdu.add_keyword(std::string("CRPIX2"),std::string("1"),empty_string);
    hdu.add_keyword(std::string("CRVAL2"),freq_start,empty_string);
    hdu.add_keyword(std::string("CDELT2"),delta_freq,empty_string);
   
    hdu.add_keyword(std::string("CTYPE1"),std::string("TIME"),empty_string);
    hdu.add_keyword(std::string("CUNIT1"),std::string("sec"),empty_string);
    hdu.add_keyword(std::string("CRPIX1"),1.00,empty_string);
    hdu.add_keyword(std::string("CRVAL1"),0.00,empty_string);
    hdu.add_keyword(std::string("CDELT1"),delta_time,empty_string);
    hdu.add_keyword(std::string("INTTIME"),delta_time,empty_string);
    
    fitsImage.add_HDU(hdu);
    fitsImage.write();
}
