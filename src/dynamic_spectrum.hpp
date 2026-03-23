#ifndef __DYNAMIC_SPECTRUM__
#define __DYNAMIC_SPECTRUM__
#include <cstring>
#include <stdexcept>
#include <images.hpp>
#include <memory_buffer.hpp>

class DynamicSpectrum : public MemoryBuffer<float> {

    private:
    size_t n_timesteps {0}, n_channels {0}, current_offset {0}, batch_size {0};
    // image pixel coordinates this dynamic spectrum refers to
    float freq_start {138.9}, delta_freq{0.04}, delta_time{0.02};
    
    public:
    int x {0}, y {0};
    
    DynamicSpectrum(size_t n_timesteps, size_t n_channels, size_t batch_size,  int x, int y) : \
            MemoryBuffer {n_timesteps * n_channels, MemoryType::MANAGED}, n_timesteps {n_timesteps}, \
            n_channels {n_channels}, x {x}, y {y}, batch_size {batch_size} {
        std::memset(this->data(), 0, sizeof(float) * this->size());
    }

    void add_images(Images& images);

    void increase_offset() {current_offset += batch_size;};

    void to_fits_file(std::string filename);
    
    // setters:
    inline void set_freq_start(float freqstart_param){ freq_start = freqstart_param; }
    inline void set_delta_freq(float delta_freq_param){ delta_freq = delta_freq_param; }
    inline void set_delta_time(float delta_time_param){ delta_time = delta_time_param; }
    
    // getters:
    inline int get_ntimesteps(){ return n_timesteps; }
    inline int get_nchannels(){ return n_channels; }
    inline int get_current_offset(){ return current_offset; }
    inline int get_batch_size(){ return batch_size; }
    inline float getfreq_start(){ return freq_start; }
    inline float get_delta_freq(){ return delta_freq; }
    inline float get_delta_time(){ return delta_time; }
};

#endif