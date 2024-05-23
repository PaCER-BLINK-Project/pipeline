#include <iostream>
#include <vector>
#include <string>
#include <unistd.h>
#include <stdexcept>
#include <cstring>
#include <sstream>
#include <memory>
#include "../src/pipeline.hpp"

#include <pacer_imager.h>
// profiling :
#include <mystring.h>
#include <myparser.h>
#include <mydate.h>
#include <mystrtable.h>
#include <pacer_imager_defs.h>
#include "../src/files.hpp"



void print_program_options(const blink::ProgramOptions& opts);
void parse_program_options(int argc, char** argv, blink::ProgramOptions& opts);
void print_help(std::string exec_name, blink::ProgramOptions& opts);



int main(int argc, char **argv){
    blink::ProgramOptions opts;

     if(argc < 2){
        print_help(argv[0], opts);
        exit(0);
    }
    try {
        parse_program_options(argc, argv, opts);
    } catch (std::invalid_argument& ex){
        std::cerr << ex.what() << std::endl;
        exit(1);
    }
     // Create output directory if exists
    if(opts.outputDir != "." && ! blink::dir_exists(opts.outputDir)){
        if(!blink::create_directory(opts.outputDir)){
            std::cerr << "Impossible to create the output directory." << std::endl;
            exit(1);
        }
    }
    print_program_options(opts);


    // TODO : replace all these options to just ProgramOptions& opts - not doing it now !
    blink::Pipeline pipeline  {
        opts.nChannelsToAvg, opts.integrationTime, opts.reorder, opts.szCalibrationSolutionsFile.length() > 0,
        opts.szCalibrationSolutionsFile, opts.ImageSize, opts.MetaDataFile,
        opts.szAntennaPositionsFile, opts.MinUV, opts.bPrintImageStatistics, opts.szWeighting,
        opts.outputDir, opts.bZenithImage, opts.FrequencyMHz, opts.FOV_degrees, opts.inputDataType, 
        opts.fUnixTime, opts.ApplyCalibrationInImager, opts.gFlaggedAntennasList, opts.outputDir
    };

    bool on_gpu = num_available_gpus() > 0;
    if(opts.input_files.size() == 1){
        std::string& filename = opts.input_files[0];
        ObservationInfo obs_info = parse_mwa_phase1_dat_file_info(filename);
        obs_info.coarse_channel_index = opts.coarseChannelIndex;
        unsigned int integration_steps {static_cast<unsigned int>(opts.integrationTime / obs_info.timeResolution)};
        auto volt = Voltages::from_dat_file(filename, obs_info, integration_steps, on_gpu);
        pipeline.run(volt ,opts.FreqChannelToImage);
    }else{
        auto observation = parse_mwa_dat_files(opts.input_files);
        for (auto& one_second_data : observation) {
            /* the std::vector memory allocator progressively allocates larger chunk of memory to
             * make space for new elements, copying existing elements to the new memory location.
             * In turn the MemoryBuffer copy contructor is called multiple times for the same object,
             * wasting time and memory. The solution is to use smart pointers to Voltages object.
             * In this way, smart pointers are the ones being copied over and over, instead of
             * the object they point to.
             */
            std::vector<std::shared_ptr<Voltages>> voltages;
            for(auto& dat_file : one_second_data){
                std::string& filename {dat_file.first};
                ObservationInfo obs_info {dat_file.second};
                obs_info.coarse_channel_index = static_cast<unsigned int>(voltages.size());
                unsigned int integration_steps {static_cast<unsigned int>(opts.integrationTime / obs_info.timeResolution)};
                auto volt = Voltages::from_dat_file(filename, obs_info, integration_steps, on_gpu);
                voltages.push_back(std::make_shared<Voltages>(std::move(volt)));
            }
            pipeline.run(voltages ,opts.FreqChannelToImage);
        }
    }
}



void print_help(std::string exec_name, blink::ProgramOptions& opts ){
    std::cout << "\n" << exec_name << " -t <int time> [-c <cnls>] [-o <outdir>] DATFILE1 [DATFILE2 [DATFILE3 [...]]]\n"
    "\nProgram options:\n-------------\n"
    "\nCorrelator options:\n"
    "\t-t <integration time>: duration of the time interval to integrate over. Accepts a timespec (see below).\n"
    "\t-c <channels to average>: number of contiguous frequency channels to average. Must be >= 1.\n"
    "\t\t Default is 1, that is, no averaging.\n"
    "\t-o <output directory>: path to a directory where to save output files. If the directory does\n"
    "\t\t not exist, it will be created. Default is current directory.\n"
    "\t-s <solutions file>: calibration solutions file. If specified, calibration solutions will be applied to the\n"
    "\t-b : apply calibration straight after correlation (not in the imager - default)\n"
    "\t\t correlated visibilities.\n"
    "\n"    
    "Time specification (timespec)\n-----------------------------\n"
    "A timespec is a convenient way of specifing a time duration. It is made up of a integer number\n"
    "followed by a unit of time; for instance '2ms' is a timespec representing 2 milliseconds.\n"
    "Valid unit of times are: 'ms', 'cs', 'ds', and 's'.\n"
    "\n\n\n\nImager options:\n"
    "\t-p POSTFIX : default is not postfix it's added to basename specified with VISIBILITY_FITS_BASENAME\n"
    "\t-f FREQ_MHz     : frequency in MHz [default 159.375 MHz]\n" // should be opts.FrequencyMHz but I am not used to cout 
    "\t-F FoV[deg]     : field of view in degrees [default 180 degree]\n" // should be opts.FOV_degrees but I am not used to cout
    "\t-w WEIGHTING    : change weighting schema N - natural, U - uniform [default N]\n" // should be opts.szWeighting    
    "\t-m MIN_UV_DISTANCE : minimum UV distance in wavelengths for a baseline to be included [default -1000]\n" // change to opts.MinUV
    "\t-a antenna_positions.txt : text file with antenna positions in a format : AntName X[m] Y[m] Z[m]\n"
    "\t-A flagged antennas list string (coma separated, e.g. 1,2,3,4) [default empty]\n"
    "\t-n IMAGE_SIZE : single value N for image size N x N pixels [default 128]\n" // change to opts.ImageSize
    "\t-Z : image phase centered at zenith, re-calculation of UVW is not required\n" 
    "\t-S : print statistics of the final sky image (mean,median,rms,rms_iqr etc)\n"
    "\t-M META_DATA_FILE : name of meta data file\n"    
    "\t-U unixtime : unix time of the start of the data [default none -> use current time]\n"
    "\t-v VERBOSITY/DEBUG level [default ??? ]\n" // TO-ADD CPacerImager::m_ImagerDebugLevel but I am not used to cout 
    "\t-V FILE_SAVE_LEVEL : to control number of FITS files saved [default ???]\n" // TO-ADD CPacerImager::m_SaveFilesLevel I am not used to cout 
    "\t-C frequency_channel to image [default -1 - means all]\n"
    "\t-G : apply geometric correction\n"
    "\t-L : apply cable correction\n"  
    "\t"
    << std::endl;
}



void parse_program_options(int argc, char** argv, blink::ProgramOptions& opts){
    // default values
    opts.outputDir = std::string {"."};
    opts.nChannelsToAvg = 1;
    opts.integrationTime = -1.0; // No default!
    opts.inputDataType = blink::DataType::MWA;
    opts.szCalibrationSolutionsFile = "";
    opts.MinUV = -1000;
    opts.bZenithImage = false;
    opts.szWeighting = "N";
    opts.ImageSize = 180;
    opts.reorder = false;
    opts.bPrintImageStatistics = false;
    opts.FreqChannelToImage = -1; // image all channels 
    opts.ApplyCalibrationInImager = true; // default to apply calibration in the imager
    
    // default debug levels :
    CPacerImager::SetFileLevel(SAVE_FILES_FINAL);
    CPacerImager::SetDebugLevel(IMAGER_WARNING_LEVEL);

    const char *options = "rt:c:o:a:M:Zi:s:f:F:n:U:v:w:V:C:GLA:b:";
    int current_opt;
    while((current_opt = getopt(argc, argv, options)) != - 1){
        switch(current_opt){
            case 'r': {
                opts.reorder = 1;
                break;
            }
            case 'a': {
               if( optarg && strlen(optarg) ){
                  opts.szAntennaPositionsFile = optarg;
               }
               break;
            }
            
            case 'A': {
               if( optarg && strlen(optarg) ){
                  opts.gFlaggedAntennasListString = optarg;
                  
               }
               break;
            }
            
            case 'b' : {               
               opts.ApplyCalibrationInImager = false;
               opts.coarseChannelIndex = atoi(optarg);
               break;
            }

            case 'C': {
               if( optarg && strlen(optarg) ){
                  opts.FreqChannelToImage = atol( optarg );
               }
               break;
            }


            case 't': {
                opts.integrationTime = parse_timespec(optarg);
                if(opts.integrationTime == 0) throw std::invalid_argument("Non positive integration time specified.");
                break;
            }
            case 'o' : {
                opts.outputDir = std::string {optarg};
                break;
            }
            case 'L' : {
                CImagerParameters::m_bApplyCableCorr = true;
                break;
            }
            case 'G' : {
                CImagerParameters::m_bApplyGeomCorr = true;
                break;
            }            
            case 'M': {
               opts.MetaDataFile = std::string {optarg};
               // imager.m_ImagerParameters.m_bConstantUVW = false; // when Meta data file is provided it assumes that it will pointed observation (not all sky)
               opts.bZenithImage = false;
               break;
            }
            case 'c': {
                opts.nChannelsToAvg = atoi(optarg);
                if(opts.nChannelsToAvg < 1) throw std::invalid_argument("Value for number of channels to average must be at least 1.");
                break;
            }
            
            case 'i': {
                if(!strcmp("eda2", optarg)) opts.inputDataType = blink::DataType::EDA2;
                else if(strcmp("mwa", optarg)){
                    std::stringstream ss;
                    ss << "Unrecognised data type: '" << optarg  << "'.";
                    throw std::invalid_argument(ss.str());
                }
                break;
            }
            
            case 'Z': {
               opts.bZenithImage = true;
               // imager.m_ImagerParameters.m_bConstantUVW = true; // constant UVW (not re-calculated as a function of time due to pointing)
               // imager.m_ImagerParameters.m_bConstantUVW = false; // when Meta data file is provided it assumes that it will pointed observation (not all sky)
               opts.bZenithImage = false;
               break;
            }
            case 's' : {
                opts.szCalibrationSolutionsFile = std::string {optarg};
                break;

            }
            case 'f' : {
                opts.FrequencyMHz = atof(optarg);
                break;
            }
            case 'F': {
                opts.FOV_degrees = atof(optarg);
                break;
            }
            case 'n': {
                opts.ImageSize = atoi(optarg);
                break;
            }
            
            case 'U' : {
               if( optarg && strlen(optarg) ){
                   opts.fUnixTime = atof( optarg );
               }
               break;
            }
            
            case 'w' : {
               if( optarg ){
                  opts.szWeighting = optarg;
               }
               break; 
            }
            
            case 'v' : {
               if( optarg ){
                   CPacerImager::SetDebugLevel( atol( optarg ) );                  
               }
               break; 
            }
            
            case 'V' : {
               if( optarg ){
                  CPacerImager::SetFileLevel( atol( optarg ) );
               }
               break; 
            }

            default : {
                std::stringstream ss;
                ss << "Unrecognised option: '" << static_cast<char>(optopt) << "'.";
                throw std::invalid_argument(ss.str());
            }
        }
    }
    for(; optind < argc; optind++) opts.input_files.push_back({argv[optind]});
    
    // options validation
    if(opts.integrationTime <= 0) throw std::invalid_argument("You must specify a value for the integration time.");
    if(opts.input_files.size() == 0) throw std::invalid_argument("No input file specified.");
    
    if( strlen(opts.gFlaggedAntennasListString.c_str()) > 0 ){
       MyParser pars=opts.gFlaggedAntennasListString.c_str();
       CMyStrTable items;
       pars.GetItemsNew( items, "," );
       for(int i=0;i<items.size();i++){
          opts.gFlaggedAntennasList.push_back( atol( items[i].c_str() ) );
       }
    }
}


void print_program_options(const blink::ProgramOptions& opts){
    std::cout << "\nRunning the correlator program with the following options:\n"
    "\t Integration time interval: " << opts.integrationTime << "s\n"
    "\t Freq. channel: " << opts.FreqChannelToImage << "\n"
    "\t Number of channels to average: " << opts.nChannelsToAvg << "\n"
    "\t Output directory: " << opts.outputDir << "\n" 
    "\t Data type: " << opts.inputDataType  << "\n"
    "\t Calibration in imager: " << opts.ApplyCalibrationInImager << "\n"
    "\t Calibration file: " << opts.szCalibrationSolutionsFile << "\n"    
    "\t Zenith image: " << opts.bZenithImage << std::endl;
}
