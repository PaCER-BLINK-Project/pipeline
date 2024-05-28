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
    double minUV, bool printImageStats, std::string szWeighting, std::string outputDir, bool bZenithImage, double frequencyMHz, double FOV_degrees, blink::DataType inputType, double fUnixTime, bool b_calibrate_in_imager, vector<int>& flagged_antennas, std::string& output_dir){

    // set imager parameters according to options :    
    // no if here - assuming always true :
    this->imager.m_ImagerParameters.m_bConstantUVW = true; 
    this->channels_to_avg = nChannelsToAvg;
    this->integration_time = integrationTime;
    this->calibrate = calibrate;
    this->calibration_solutions_file = solutions_file;
    this->calibrate_in_imager = b_calibrate_in_imager;
    this->imageSize = imageSize;
    this->MetaDataFile = metadataFile;
    this->MinUV = minUV;
    this->bPrintImageStatistics = printImageStats;
    this->szWeighting = szWeighting;
    this->bZenithImage = bZenithImage;
    this->inputDataType = inputType;
    this->frequencyMHz = frequencyMHz;
    this->FOV_degrees = FOV_degrees;
    this->reorder = reorder;
    this->output_dir = output_dir;
    imager.m_ImagerParameters.SetGlobalParameters(szAntennaPositionsFile.c_str(), bZenithImage); // Constant UVW when zenith image (-Z)
    imager.m_ImagerParameters.m_szOutputDirectory = outputDir.c_str();
    imager.m_ImagerParameters.m_fCenterFrequencyMHz = frequencyMHz;


    if(strlen(metadataFile.c_str())){
        imager.m_ImagerParameters.m_MetaDataFile = metadataFile.c_str();
        imager.m_ImagerParameters.m_bConstantUVW = false; // when Meta data file is provided it assumes that it will pointed observation (not all sky)
        this->bZenithImage = false;
    }
    
    imager.m_ImagerParameters.m_fUnixTime = fUnixTime;
    if( this->calibrate_in_imager ){
       printf("DEBUG : setting calibration solution file in imager (to apply calibration there)\n");
       imager.m_CalibrationSolutions.m_filename = solutions_file.c_str();
    }else{
       printf("DEBUG : setting calibration solution file in imager to empty value (apply cal. in blink_pipeline)\n");
       imager.m_CalibrationSolutions.m_filename = "";
    }


    // // overwrite with number of antennas according to the list :
    // if( imager.m_MetaData.m_AntennaPositions.size() > 0 ){
        // Cristian's comment: this can be dangerous for the correlator
    //     opts.obsInfo.nAntennas = imager.m_MetaData.m_AntennaPositions.size();
    // }

    // // set antenna flags (if there are):
    // if( opts.szFlaggedAntennasList.size() > 0 ){
    //     imager.SetFlaggedAntennas( opts.szFlaggedAntennasList );
    // }
   // imager.SetFileLevel(SAVE_FILES_FINAL);
   // imager.SetDebugLevel(IMAGER_WARNING_LEVEL);

   // TODO : de-hardcode , add option for a list of flagged antennas !
   // vector<int> flagged_antennas;
   // flagged_antennas.push_back(80);
   // flagged_antennas.push_back(119);
   // imager.SetFlaggedAntennas( flagged_antennas );
   imager.Initialise();
   
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
   
   auto xcorr = cross_correlation(input, channels_to_avg);
   
   if(reorder){
      auto mapping = get_visibilities_mapping(this->MetaDataFile);
	   xcorr = reorder_visibilities(xcorr, mapping);
   }

   if( calibrate_in_imager ){ 
      printf("DEBUG : calibration will be applied in imager (not in blink::Pipeline::run)\n");
   }else{
      if( calibrate ){ // disabled for now before I check other things
         std::cout << "Calibration is being applied in the pipeline ( coarse channel index = " << obsInfo.coarse_channel_index << ")." << std::endl;
         auto sol = CalibrationSolutions::from_file(this->calibration_solutions_file);
         apply_solutions(xcorr, sol, obsInfo.coarse_channel_index);
      }
   }

   // TODO : keep the loop so that it change be parallelised using OpenMP
   // xcorr.nFrequencies = 1;
   printf("DEBUG : imaging %d intervals and %d frequency channels\n",int(xcorr.integration_intervals()),xcorr.nFrequencies);
   for(int integrationInterval {0}; integrationInterval < xcorr.integration_intervals(); integrationInterval++){
      for(int frequency {0}; frequency < xcorr.nFrequencies; frequency++){
         if ( frequency == freq_channel || freq_channel < 0 ){
            // if freq_channel is specified (processing of a single channel) -> use frequencyMHz as the center frequency of this channel
            // otherwise, calculate frequency assuming frequencyMHz is the start frequency of the entire band:
            // I am including cotter-compatibility option here to be able to compare with standard SMART pipeline processing:
            double coarse_channel_central_freq_MHz = obsInfo.coarseChannel * 1.28;
            double channel_frequency_MHz = frequencyMHz;
            if( freq_channel < 0 ){ // processing of all channels requires re-calculation of frequency :
               bool cotter_compatible=true;
               double fine_ch_bw = 0.04, coarse_ch_bw=1.28;
               if( cotter_compatible ){
                  // see for example awk commands in experiments/image_mwa_obsid1276619416_allch.sh 
                  // frequencyMHz is assumed to be center frequency of coarse channel - 0.64 -> lower edge of the MWA coarse channel:
                  // TODO : get this information from xcorr structure or metedata, otherwise it will not work for EDA2 etc:                  
                  //        it looks to me that the required fields need to be added there first
                  channel_frequency_MHz = coarse_channel_central_freq_MHz - coarse_ch_bw/2.00 + fine_ch_bw*frequency; // cotter has a bug ? - does not add half of fine channel to calculate channel frequency 
               }else{
                  channel_frequency_MHz = coarse_channel_central_freq_MHz - coarse_ch_bw/2.00 + fine_ch_bw*frequency + fine_ch_bw/2.00;
               }          
               imager.m_ImagerParameters.m_fCenterFrequencyMHz = channel_frequency_MHz; // update parameter in the imager too to make sure frequnecies are consistent 
               
               char szOutDir[64];
               sprintf(szOutDir,"%s/%ld/%d/%03d", output_dir.c_str(), obsInfo.startTime, obsInfo.coarseChannel, frequency);
               imager.m_ImagerParameters.m_szOutputDirectory = szOutDir;
               if(calibrate_in_imager){
                  char szChannelSolutionFile[1024];
                  // use calibration_solutions_file.c_str() as a basename and add channel 
                  sprintf(szChannelSolutionFile,"%s_chan%03d_xx.txt",calibration_solutions_file.c_str(),(obsInfo.coarse_channel_index * xcorr.nFrequencies + frequency)); // WARNING : this is test version hence XX hardcoded (needs to be for both XX and YY)
                  std::cout << "Using solution file " << szChannelSolutionFile << std::endl;
                  imager.UpdateCalSolFile( szChannelSolutionFile );
               }         
            }
         
            //----------------------------------------------------------------
            // MS : 20230914 - temporary test code:
            // CBgFits vis_re(128,128),vis_im(128,128);
            // char szOutPutFits[1024];
            // ConvertXCorr2Fits( xcorr, vis_re, vis_im, integrationInterval, frequency, imager.m_ImagerParameters.m_szOutputDirectory.c_str() );
            // sprintf(szOutPutFits,"%s/re.fits",imager.m_ImagerParameters.m_szOutputDirectory.c_str());
            // vis_re.WriteFits( szOutPutFits );
            // sprintf(szOutPutFits,"%s/im.fits",imager.m_ImagerParameters.m_szOutputDirectory.c_str());
            // vis_im.WriteFits( szOutPutFits );
            // //----- end of temporary test code

            printf("DEBUG : starting imager using xcorr structure ( frequency = %d , frequencyMHz = %.6f [MHz] -> channel frequency = %.6f [MHz] )\n",frequency,frequencyMHz,channel_frequency_MHz);
            char szOutImage[64];
            sprintf(szOutImage,"test_image_time%06d_ch%05d",integrationInterval,frequency);
            imager.run_imager( xcorr, integrationInterval, frequency, channel_frequency_MHz, imageSize, FOV_degrees, 
                                 MinUV, true, true, szWeighting.c_str(), szOutImage, false);
         }
      }
   }
}
