#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "peak_finding.hpp"

struct Candidate {
    int x;
    int y;
    float dm;
    int time_step;
    float value;
    float mean;
    float stdev;
};

struct TimeSeriesCoord {
    int x;
    int y;
    int dm_idx;
    float dm;
};

void peakfinding_simple_avg(float *dm_starttime, int side_size, std::vector<float>& dm_list, int window_start_idx, int buffer_size, 
        int table_size, int global_offset, const std::vector<float>& norm_factors, float SNR, std::string filename){
    
    int n_dms = dm_list.size();
    unsigned long long  width = n_dms * table_size * side_size;
    unsigned long long  cell_size = n_dms * table_size;

    std::vector<Candidate> candidates;
    std::vector<TimeSeriesCoord> time_series;

    #pragma omp parallel for collapse(2)
    for(int x = 0; x < side_size; x++){
        for(int y = 0; y < side_size; y++){
            for(int dm {0}; dm < n_dms; dm++){
                float avg = 0.0;
                float stdev = 0.0;
                bool save_series {false};
                for(int i {0}; i < buffer_size; i++){
                    int idx = (window_start_idx + i) % table_size;
                    dm_starttime[y * width + x*cell_size + dm * table_size + idx] /= norm_factors[dm];
                    avg += dm_starttime[y * width + x*cell_size + dm * table_size + idx];
                }
                avg /= buffer_size;
                for(int i {0}; i < buffer_size; i++){
                    int idx = (window_start_idx + i) % table_size;
                    float val =  dm_starttime[y * width + x*cell_size + dm * table_size + idx];
                    stdev += std::pow(val - avg, 2);
                }
                stdev = std::sqrt(stdev / buffer_size);
                for(int i {0}; i < buffer_size; i++){
                    int idx = (window_start_idx + i) % table_size;
                    float val =  dm_starttime[y * width + x*cell_size + dm * table_size + idx];
                    if(((val - avg) / stdev) >= SNR){
                        #pragma omp critical 
                        {   
                            Candidate c {
                                .x = x,
                                .y = y,
                                .dm = dm_list[dm],
                                .time_step = i,
                                .value = val,
                                .mean = avg,
                                .stdev = stdev
                            };
                            candidates.push_back(c);
                            if(!save_series){
                                save_series = true;
                                time_series.push_back({x, y, dm, dm_list[dm]});
                            }
                        }
                    }
                } 
            }
        }
    }
    std::ofstream out_file(filename, std::ios::app | std::ios::binary);
    #define write_to_file(X) out_file.write(reinterpret_cast<char*>(&X), sizeof(X))

    write_to_file(global_offset);
    write_to_file(side_size);
    write_to_file(n_dms);
    write_to_file(buffer_size);
    int n_candidates = static_cast<int>(candidates.size());
    int n_series = static_cast<int>(time_series.size());
    write_to_file(n_candidates);
    write_to_file(n_series);
    for(int i {0}; i < n_candidates; i++)
        write_to_file(candidates[i]);
    for(int j {0}; j < n_series; j++){
        write_to_file(time_series[j].x);
        write_to_file(time_series[j].y);
        write_to_file(time_series[j].dm);
        for(int i {0}; i < buffer_size; i++){
            int idx = (window_start_idx + i) % table_size;
            float val =  dm_starttime[time_series[j].y * width + time_series[j].x * cell_size + \
                time_series[j].dm_idx * table_size + idx];
            write_to_file(val);
        }
    }
    out_file.close();
}
