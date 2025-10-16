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
#include "gpu/dedispersion_gpu.hpp"
#include "peak_finding.hpp"
#include <gpu_macros.hpp>
#include "rfi_flagging.hpp"
#include "gpu/rfi_flagging_gpu.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif
// TODO: we could use an observation info structure here to pass values
blink::Pipeline::Pipeline(unsigned int nChannelsToAvg, double integrationTime, bool reorder, bool calibrate,
    std::string solutions_file, int imageSize, std::string metadataFile, float oversampling_factor,
    std::string szAntennaPositionsFile, double minUV, bool printImageStats, std::string szWeighting, 
    std::string outputDir, bool bZenithImage, double FOV_degrees, bool averageImages, Polarization pol_to_image,
    vector<int>& flagged_antennas, bool change_phase_centre, double ra_deg, double dec_deg, Dedispersion& dedisp_engine,
    float rfi_flagging, std::string& output_dir, std::string& postfix) : dedisp_engine {dedisp_engine}, total_power(integrationTime*1000.00) {

    gpuGetDeviceCount(&num_gpus);
    if(num_gpus == 0){
        std::cerr << "No GPU device detected. Aborting.." << std::endl;
        throw std::exception();
    }
    imager.resize(num_gpus);
    mapping.resize(num_gpus);
    cal_sol.resize(num_gpus);

    // set imager parameters according to options :    
    // no if here - assuming always true :
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
    this->postfix = postfix;
    this->bZenithImage = false;
    this->rfi_flagging = rfi_flagging;

    if(calibrate) cal_sol[0] = CalibrationSolutions::from_file(this->calibration_solutions_file);
    if(reorder) mapping[0] = get_visibilities_mapping(this->MetaDataFile);
    for(int i {1}; i < num_gpus; i++){
        if(calibrate) cal_sol[i] = cal_sol[0];
        if(reorder) mapping[i] = mapping[0];
    }
    for(int i {0}; i < num_gpus; i++){
        gpuSetDevice(i);
        imager[i] = new CPacerImagerHip {metadataFile, imageSize, flagged_antennas, averageImages,
            pol_to_image, oversampling_factor, MinUV, szWeighting.c_str()};
        if(change_phase_centre) imager[i]->m_MetaData.set_phase_centre(ra_deg, dec_deg);
        mapping[i].to_gpu();
        cal_sol[i].to_gpu();
    }
}


void blink::Pipeline::run(const std::vector<std::shared_ptr<Voltages>>& inputs){
    if(dedisp_engine.is_initialised() && dedisp_engine.buffer_is_full())
        dedisp_engine.process_buffer();
   
   #pragma omp parallel for num_threads(num_gpus)
   for(int i = 0; i < inputs.size(); i++){
        const auto& input = inputs[i];
        #ifdef _OPENMP
        run(*input, omp_get_thread_num());
        #else
        run(*input, 0);
        #endif
   }
   if(dedisp_engine.is_initialised()) dedisp_engine.increase_offset();
   if(DynamicSpectra.size() > 0 ){
      for( auto ds : DynamicSpectra ){
         ds->increase_offset();
      }
   }
}

void blink::Pipeline::run(const Voltages& input, int gpu_id){

    gpuSetDevice(gpu_id);
    const ObservationInfo& obsInfo {input.obsInfo};
    unsigned int nIntegrationSteps {static_cast<unsigned int>(integration_time / obsInfo.timeResolution)};

    std::cout << "Correlating voltages (OBSID = " << obsInfo.id << ", Coarse Channel = " << obsInfo.coarseChannel << ") .." << std::endl;

    high_resolution_clock::time_point corr_start = high_resolution_clock::now();
    auto xcorr = cross_correlation(input, channels_to_avg);
    high_resolution_clock::time_point corr_end = high_resolution_clock::now();
    duration<double> corr_dur = duration_cast<duration<double>>(corr_end - corr_start);
    std::cout << "Cross correlation took " << corr_dur.count() << " seconds." << std::endl;

    if(reorder){
        xcorr = reorder_visibilities(xcorr, mapping[gpu_id]);
    }

    if( calibrate ){
        std::cout << "Calibration is being applied in the pipeline ( coarse channel index = " << obsInfo.coarse_channel_index << ")." << std::endl;
        apply_solutions(xcorr, cal_sol[gpu_id], obsInfo.coarse_channel_index);
    }

    std::cout << "Running imager.." << std::endl;
    auto images = imager[gpu_id]->run(xcorr);
   
    if(rfi_flagging > 0){
        // images.to_cpu();
        std::cout << "Applying RFI flagging..." << std::endl;
        //flag_rfi_cpu(images, rfi_flagging, 50, false);
        flag_rfi_gpu(images, rfi_flagging);
        clear_flagged_images_gpu(images);
    }
    if(DynamicSpectra.size() > 0 ) {
        std::cout << "Adding images to dynamic spectrum.." << std::endl;
        images.to_cpu();
        // TODO add GPU implementation
        for( auto ds : DynamicSpectra ){
           ds->add_images(images);
        }
    }

    if(false){
       // TODO - set offset and ntimes and start_time to make sure the new portion of data is processed
       int offset=0; // use DynamicSpectra.current_offset use this one 
       int ntimes=1; // use DynamicSpectra.batch_size ??? - Cristian : how many timesteps are processed each time (20ms) - 50 in our case !!!
       int start_time=0;
       total_power.calc( *(DynamicSpectra[0]), true, true, offset, ntimes, start_time );
       // determine flags for images in the same time step
       // process(DynamicSpectra[0], images);

       // clear_flagged_images_gpu(images);
    }

    if(DynamicSpectra.size() == 0 && !dedisp_engine.is_initialised()){
        std::cout << "Saving images to disk..." << std::endl;
        images.to_fits_files(output_dir);
    }        

    if(dedisp_engine.is_initialised()){
        int top_freq_idx = (obsInfo.coarse_channel_index + 1) * images.n_channels - 1;
        high_resolution_clock::time_point dedisp_start = high_resolution_clock::now();
        dedisp_engine.compute_partial_dedispersion_gpu(images, top_freq_idx);
        high_resolution_clock::time_point dedisp_end = high_resolution_clock::now();
        duration<double> dedisp_dur = duration_cast<duration<double>>(dedisp_end-dedisp_start);
        std::cout << "Dedispersion took " << dedisp_dur.count() << " seconds." << std::endl;
    }
}

void blink::Pipeline::add_dynamic_spectrum(std::shared_ptr<DynamicSpectrum> p)
{
   DynamicSpectra.push_back(p);
}

void blink::Pipeline::save_dynamic_spectra() 
{
   for( auto ds : DynamicSpectra ){
      char filename[2048];
      sprintf(filename,"%s/dynamic_spectrum_%05d_%05d.fits",output_dir.c_str(),ds->x,ds->y);
      
      ds->to_fits_file(filename);
      
      cout << "Saved dynamic spectrum " << filename << endl;      
   }
}
