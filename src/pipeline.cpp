#include <iostream>
#include <vector>
#include <string>
#include <unistd.h>
#include <stdexcept>
#include <cstring>
#include <sstream>
#include <memory>
#include <chrono>
// #include <filesystem>

using namespace std::chrono;

#include <correlation.hpp>
#include <calibration.hpp>
#include <mapping.hpp> // to apply cotter mapping
#include <pacer_imager_defs.h>
#include "pipeline.hpp"
#include "files.hpp"

// MS : 20230914 - temporary test code:
// TODO : remove this function in the future use complex<double> 
void ConvertXCorr2Fits(Visibilities& xcorr, CBgFits& vis_re, CBgFits& vis_im, int time_step, int fine_channel, const char* szOutputDir="./" )
{

   int n_ant = xcorr.obsInfo.nAntennas;
   int n_corrs = 4; // 4 correlation products : XX XY YX YY 
   int n_baselines = n_ant*(n_ant+1)/2;
   ObservationInfo& obsInfo = xcorr.obsInfo;

   const size_t matrixSize = n_baselines * obsInfo.nPolarizations * obsInfo.nPolarizations;
   const size_t nIntervals  = (obsInfo.nTimesteps); // TODO  + voltages.nIntegrationSteps - 1) / voltages.nIntegrationSteps;
   unsigned int nChannelsToAvg = 1; // TODO : verify
   const size_t nOutFrequencies = obsInfo.nFrequencies / nChannelsToAvg;
   const size_t nValuesInTimeInterval = matrixSize * nOutFrequencies;
   const size_t outSize = nValuesInTimeInterval * nIntervals;
   unsigned int avgCh = fine_channel / nChannelsToAvg;

   

   // size_t outIndex {interval * nValuesInTimeInterval + avgCh * matrixSize + baseline * obsInfo.nPolarizations * obsInfo.nPolarizations
   //                             + p1*obsInfo.nPolarizations + p2};
   // TODO !!!
   // assuming single timestep for now :  
   // assuming ordering of data as Cristian told me during the meeting :
   // Ant11        | Ant 12       | Ant 13       | ...
   // XX XY YX YY  | XX XY YX YY  | XX XY YX YY  | ...
   int index = 0;
 
   vis_re.SetNaN();
   vis_im.SetNaN();
   // using at :
   // Complex<float> *at_float(Visibilities& vis, unsigned int interval, unsigned int frequency, unsigned int a1, unsigned a2)
   for(int i=0;i<n_ant;i++){ // loop over ant1          
     for(int j=0;j<=i;j++){ // loop over ant2 
        // auto& vis = xcorr.data[idx];
//        std::complex<double>* vis = at( xcorr, time_step, fine_channel, i, j );
        std::complex<float>* vis = xcorr.at( time_step, fine_channel, i, j );
       
        // current signs of imaginary work ok with Cristian's re-mapping:
        vis_re.setXY(j,i,float(vis[0].real()));
        vis_im.setXY(j,i,(+1)*float(vis[0].imag())); // was - 
        vis_re.setXY(i,j,float(vis[0].real()));
        vis_im.setXY(i,j,(-1)*float(vis[0].imag())); // was +
     }     

     index += (n_ant-i);
  }

  char szOutPutFits[1024];
  sprintf(szOutPutFits,"%s/test_vis_re.fits",szOutputDir);  
  vis_re.WriteFits( szOutPutFits );
  sprintf(szOutPutFits,"%s/test_vis_im.fits",szOutputDir);
  vis_im.WriteFits( szOutPutFits );
}
//-----------------------------------------------------------------------------------------------


blink::Pipeline::Pipeline(unsigned int nChannelsToAvg, double integrationTime, bool reorder, bool calibrate, std::string solutions_file, int imageSize, std::string metadataFile, std::string szAntennaPositionsFile,
    double minUV, bool printImageStats, std::string szWeighting, std::string outputDir, bool bZenithImage, double FOV_degrees, blink::DataType inputType, bool averageImages, vector<int>& flagged_antennas, std::string& output_dir){

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
    this->inputDataType = inputType;
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
   
 /* The following must be done in the imager
   if( opts.bChangePhaseCentre ){
      double obsid = -1;
      if( CImagerParameters::m_bAutoFixMetaData ){
         obsid = CObsMetadata::ux2gps( imager.m_ImagerParameters.m_fUnixTime );
      }

      imager.m_MetaData.set_radec( obsid, opts.fRAdeg, opts.fDECdeg );
      printf("DEBUG : set RADEC of phase centre to (%.8f,%.8f) at obsid = %.2f\n",opts.fRAdeg,opts.fDECdeg,obsid);
   }
   */
   
   // setting flagged antennas must be called / done after reading METAFITS file:
   if( flagged_antennas.size() > 0 ){
       imager.SetFlaggedAntennas( flagged_antennas );
   }
}

void blink::Pipeline::save_image(const Voltages& r) {
   high_resolution_clock::time_point save_image_start = high_resolution_clock::now();
   auto images = imager.final_image(r, imageSize);
   images.to_fits_files(output_dir);
   high_resolution_clock::time_point save_image_end = high_resolution_clock::now();
   duration<double> save_image_dur = duration_cast<duration<double>>(save_image_end - save_image_start);
   std::cout << "Saving image took " << save_image_dur.count() << " seconds." << std::endl;
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
   
   auto xcorr = cross_correlation(input, channels_to_avg);
   
   if(reorder){
      auto mapping = get_visibilities_mapping(this->MetaDataFile);
	   xcorr = reorder_visibilities(xcorr, mapping);
   }

   if( calibrate ){ // disabled for now before I check other things
      std::cout << "Calibration is being applied in the pipeline ( coarse channel index = " << obsInfo.coarse_channel_index << ")." << std::endl;
      auto sol = CalibrationSolutions::from_file(this->calibration_solutions_file);
      apply_solutions(xcorr, sol, obsInfo.coarse_channel_index);
   }

   std::cout << "Running imager.." << std::endl;
   auto images = imager.run_imager(xcorr, -1, -1, imageSize, FOV_degrees, 
      MinUV, true, true, szWeighting.c_str(), output_dir.c_str(), false);
   // std::cout << "Saving images to disk..." << std::endl;
   // images.to_fits_files(output_dir);
}

