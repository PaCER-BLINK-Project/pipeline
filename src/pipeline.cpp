#include <iostream>
#include <vector>
#include <string>
#include <unistd.h>
#include <stdexcept>
#include <cstring>
#include <sstream>
#include <memory>
// #include <filesystem>

#include <correlation.hpp>
#include <calibration.hpp>
#include <mapping.hpp> // to apply cotter mapping
#include <pacer_imager_defs.h>
#include "pipeline.hpp"
#include "files.hpp"
#include "dedispersion.hpp"
#include <gpu_macros.hpp>

void blink::Pipeline::set_frequencies(const std::vector<float>& frequencies){
   this->frequencies = frequencies;
   this->delay_table = compute_delay_table(frequencies, dm_list, integration_time);
   this->sweep_size = this->delay_table[(dm_list.size() - 1) * frequencies.size()] + 1;
   this->dedisp_norm_factor = compute_normalisation_factor(delay_table);
   this->table_size = sweep_size + buffer_size;
   size_t dm_starttime_size {dm_list.size() * table_size * imageSize * imageSize};
   std::cout << "Allocating " << (dm_starttime_size * sizeof(float) / (1024.0f*1024.0f*1024.0f)) << " GiB of memory for DMARRIVAL" << std::endl;
   this->dm_starttime.allocate(dm_starttime_size);
   std::memset(dm_starttime.data(), 0, sizeof(float) * dm_starttime_size);
   std::cout << "sweep_size = " << sweep_size << ", table_size = " << table_size << std::endl;
}


// TODO: we could use an observation info structure here to pass values
blink::Pipeline::Pipeline(unsigned int nChannelsToAvg, double integrationTime, bool reorder, bool calibrate, std::string solutions_file, int imageSize, std::string metadataFile, std::string szAntennaPositionsFile,
    double minUV, bool printImageStats, std::string szWeighting, std::string outputDir, bool bZenithImage, double FOV_degrees, bool averageImages, vector<int>& flagged_antennas,  std::vector<float>& dm_list, std::string& output_dir){

    // set imager parameters according to options :    
    // no if here - assuming always true :
    this->imager.m_ImagerParameters.m_bConstantUVW = true; 
    this->channels_to_avg = nChannelsToAvg;
    this->integration_time = integrationTime;
    this->calibrate = calibrate;
    this->calibration_solutions_file = solutions_file;
    this->imageSize = imageSize;
    this->MetaDataFile = metadataFile;
    this->MinUV = minUV;
    this->bPrintImageStatistics = printImageStats;
    this->szWeighting = szWeighting;
    this->bZenithImage = bZenithImage;
    this->FOV_degrees = FOV_degrees;
    this->reorder = reorder;
    this->output_dir = output_dir;
    this->dm_list = dm_list;
    this->batch_size = static_cast<int>(1.0 / integration_time); // TODO do better than hard conding
    this->buffer_size = 50; // 2 * this->batch_size;
    imager.m_ImagerParameters.SetGlobalParameters(szAntennaPositionsFile.c_str(), bZenithImage); // Constant UVW when zenith image (-Z)
    imager.m_ImagerParameters.m_szOutputDirectory = outputDir.c_str();
    imager.m_ImagerParameters.averageImages = averageImages;


    if(strlen(metadataFile.c_str())){
        imager.m_ImagerParameters.m_MetaDataFile = metadataFile.c_str();
        imager.m_ImagerParameters.m_bConstantUVW = false; // when Meta data file is provided it assumes that it will pointed observation (not all sky)
        this->bZenithImage = false;
    }
    
   imager.Initialise(0);
   
   
   // setting flagged antennas must be called / done after reading METAFITS file:
   if( flagged_antennas.size() > 0 ){
       imager.SetFlaggedAntennas( flagged_antennas );
   }
}


void blink::Pipeline::run(const std::vector<std::shared_ptr<Voltages>>& inputs){
   if(window_offset + batch_size > table_size){
      int move_ahead = sweep_size < buffer_size ? (window_offset - sweep_size) : buffer_size;
      // Time to use and clear buffer
      get_elements(dm_starttime.data(), imageSize, dm_list.size(), 0, window_start_idx, move_ahead, table_size, dedisp_norm_factor, 59, 629);
      // step 3, clear and rotate the buffer
      clear_buffer(dm_starttime.data(), imageSize, dm_list.size(), window_start_idx, table_size, move_ahead);
      window_offset -= move_ahead;
      window_start_idx = (window_start_idx + move_ahead) % table_size;
   }
   for(const auto& input : inputs){
      run(*input);
   }
   window_offset += batch_size;
}

void blink::Pipeline::run(const Voltages& input){
   const ObservationInfo& obsInfo {input.obsInfo};
   unsigned int nIntegrationSteps {static_cast<unsigned int>(integration_time / obsInfo.timeResolution)};

   std::cout << "Correlating voltages (OBSID = " << obsInfo.id << ", Coarse Channel = " << obsInfo.coarseChannel << ") .." << std::endl;
   
   high_resolution_clock::time_point corr_start = high_resolution_clock::now();
   auto xcorr = cross_correlation(input, channels_to_avg);
   high_resolution_clock::time_point corr_end = high_resolution_clock::now();
   duration<double> corr_dur = duration_cast<duration<double>>(corr_end - corr_start);
   std::cout << "Cross correlation took " << corr_dur.count() << " seconds." << std::endl;

   if(reorder){
      auto mapping = get_visibilities_mapping(this->MetaDataFile);
	   xcorr = reorder_visibilities(xcorr, mapping);
   }

   if( calibrate ){
      std::cout << "Calibration is being applied in the pipeline ( coarse channel index = " << obsInfo.coarse_channel_index << ")." << std::endl;
      auto sol = CalibrationSolutions::from_file(this->calibration_solutions_file);
      apply_solutions(xcorr, sol, obsInfo.coarse_channel_index);
   }

   std::cout << "Running imager.." << std::endl;
   auto images = imager.run_imager(xcorr, -1, -1, imageSize, FOV_degrees, 
      MinUV, true, true, szWeighting.c_str(), output_dir.c_str(), false);
   
   // std::cout << "Saving images to disk..." << std::endl;
   high_resolution_clock::time_point save_image_start = high_resolution_clock::now();
   images.to_cpu();
   
   // images.to_fits_files(output_dir);
   high_resolution_clock::time_point save_image_end = high_resolution_clock::now();
   duration<double> save_image_dur = duration_cast<duration<double>>(save_image_end - save_image_start);
   std::cout << "Copying images to CPU took " << save_image_dur.count() << " seconds." << std::endl;

   int top_freq_idx = (obsInfo.coarse_channel_index + 1) * images.nFrequencies - 1;
   high_resolution_clock::time_point dedisp_start = high_resolution_clock::now();
   compute_partial_dedispersion(images, top_freq_idx, images.nFrequencies,
      frequencies.size(), batch_size, delay_table.data(), dm_starttime.data(), dm_list.size(), table_size, window_start_idx, window_offset);
   high_resolution_clock::time_point dedisp_end = high_resolution_clock::now();
   duration<double> dedisp_dur = duration_cast<duration<double>>(dedisp_end-dedisp_start);
   std::cout << "Dedispersion took " << dedisp_dur.count() << " seconds." << std::endl;

}

