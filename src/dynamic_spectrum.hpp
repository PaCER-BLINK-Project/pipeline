#ifndef __DYNAMIC_SPECTRUM__
#define __DYNAMIC_SPECTRUM__
#include <cstring>
#include <images.hpp>
#include <memory_buffer.hpp>

class DynamicSpectrum : MemoryBuffer<float> {

    private:
    size_t n_timesteps {0}, n_frequencies {0}, current_offset {0}, batch_size {0};
    // image pixel coordinates this dynamic spectrum refers to
    
    public:
    int x {0}, y {0};
    
    DynamicSpectrum(size_t n_timesteps, size_t n_frequencies, size_t batch_size,  int x, int y) : \
            MemoryBuffer {n_timesteps * n_frequencies}, n_timesteps {n_timesteps}, \
            n_frequencies {n_frequencies}, x {x}, y {y}, batch_size {batch_size} {
        std::memset(this->data(), 0, sizeof(float) * this->size());
    }

    void add_images(const Images& images);

    void increase_offset() {current_offset += batch_size;};

    void to_fits_file(std::string filename);
};

#endif