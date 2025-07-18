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
#include <gpu_macros.hpp>



blink::Pipeline::Pipeline(unsigned int nChannelsToAvg, double integrationTime, bool reorder, bool calibrate, std::string solutions_file, int imageSize, std::string metadataFile, std::string szAntennaPositionsFile,
    double minUV, bool printImageStats, std::string szWeighting, std::string outputDir, bool bZenithImage, double FOV_degrees, bool averageImages, vector<int>& flagged_antennas, std::string& output_dir){

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


void blink::Pipeline::run(const std::vector<std::shared_ptr<Voltages>>& inputs, int freq_channel /*=-1*/ ){
   for(const auto& input : inputs){
      run(*input, freq_channel);
   }
}

void blink::Pipeline::run(const Voltages& input, int freq_channel){
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
   std::cout << "Saving images to disk..." << std::endl;
   high_resolution_clock::time_point save_image_start = high_resolution_clock::now();
   images.to_fits_files(output_dir);
   high_resolution_clock::time_point save_image_end = high_resolution_clock::now();
   duration<double> save_image_dur = duration_cast<duration<double>>(save_image_end - save_image_start);
   std::cout << "Saving image took " << save_image_dur.count() << " seconds." << std::endl;

}

