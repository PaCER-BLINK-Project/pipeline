#include <iostream>
#include <vector>
#include <string>
#include <unistd.h>
#include <stdexcept>
#include <cstring>
#include <sstream>
#include <memory>
#include <chrono>
#include <regex>
#include <array>

#include "../src/pipeline.hpp"

#include <pacer_imager.h>
// profiling :
#include <pacer_imager_defs.h>
#include <dedispersion.hpp>
#include "../src/files.hpp"
#include "../src/dynamic_spectrum.hpp"


namespace {
    std::vector<std::string> tokenize_string(std::string str, char delimiter=','){    
        std::stringstream ss(str);
        std::string token;
        std::vector<std::string> tokens;

        while (std::getline(ss, token, delimiter)) {
            tokens.push_back(token);
        }
        return tokens;
    }
}

struct ProgramOptions {
    std::string input_directory;
    int seconds_offset;
    int seconds_count;
    std::string outputDir;
    unsigned int nChannelsToAvg;
    double integrationTime;
    bool reorder;
    int coarseChannelIndex;
    // imager options:
    ObservationInfo obsInfo;
    std::string szAntennaPositionsFile;
    std::string szFlaggedAntennasListString;
    std::string szCalibrationSolutionsFile;
    std::vector<int> szFlaggedAntennasList;
    double FOV_degrees;
    std::string MetaDataFile;
    std::string ImagerOutFilePostfix;
    int ImageSize;
    double MinUV;
    bool bPrintImageStatistics;
    std::string szWeighting;
    bool bZenithImage;
    bool averageImages;
    // all-sky EDA2-like image
    int FreqChannelToImage;
    
    // change phase centre
    bool bChangePhaseCentre;
    double ra_dec[2] {0.0, 0.0};
    float SNR;
    Polarization pol_to_image;
    float oversampling_factor;
    // flagging antennas:
    string gFlaggedAntennasListString;
    vector<int> gFlaggedAntennasList;
    std::string postfix;
    string dm_list_string;
    vector<float> dm_list;
    float rfi_block_threshold;
    float rfi_history_threshold;
    vector<array<int,2>> ds_pixels;
    int n_antennas;
    int cands_per_batch;

};


void print_program_options(const ProgramOptions& opts);
void parse_program_options(int argc, char** argv, ProgramOptions& opts);
void print_help(std::string exec_name, ProgramOptions& opts);
void parse_radec_string(const std::string& radec_string, double ra_dec[]);
std::vector<float> get_frequencies(const std::vector<DatFile>& one_second, int chnls_to_avg);

int main(int argc, char **argv){
    ProgramOptions opts;

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
    int n_timesteps {static_cast<int>(1.0 / opts.integrationTime)}; // TODO: generalise to number of time steps
    blink::dedispersion::Dedispersion dedisp_engine {opts.dm_list, opts.ImageSize,
        n_timesteps, 2 * n_timesteps, opts.SNR, opts.cands_per_batch, opts.outputDir, opts.postfix};

    // TODO : replace all these options to just ProgramOptions& opts - not doing it now !
    blink::Pipeline pipeline  {
        opts.nChannelsToAvg, opts.integrationTime, opts.reorder, opts.szCalibrationSolutionsFile.length() > 0,
        opts.szCalibrationSolutionsFile, opts.ImageSize, opts.MetaDataFile, opts.oversampling_factor,
        opts.szAntennaPositionsFile, opts.MinUV, opts.bPrintImageStatistics, opts.szWeighting,
        opts.outputDir, opts.bZenithImage, opts.FOV_degrees, opts.averageImages, opts.pol_to_image,
        opts.gFlaggedAntennasList, opts.bChangePhaseCentre, opts.ra_dec[0], opts.ra_dec[1], dedisp_engine,
        opts.rfi_block_threshold, opts.rfi_history_threshold, opts.outputDir, opts.postfix
    };
   
    auto observation = parse_mwa_dat_files(opts.input_directory, opts.seconds_offset, opts.seconds_count);
    auto frequencies_list = get_frequencies(observation[0], opts.nChannelsToAvg);
    std::cout << "Will process " << observation.size() << " seconds of observation." << std::endl;

    bool dynamic_spectrum_mode_enabled = (opts.ds_pixels.size()>0);
    if(dynamic_spectrum_mode_enabled){
        for( auto ds_pixel : opts.ds_pixels ){
           auto pDynamicSpectrum = std::make_shared<DynamicSpectrum>(observation.size() * n_timesteps,
               frequencies_list.size() - 1, n_timesteps, ds_pixel[0], ds_pixel[1]);
           if( pipeline.has_dynamic_spectrum( ds_pixel[0], ds_pixel[1] ) ){
              std::cout << "Duplicate pixel " << ds_pixel[0] << "," << ds_pixel[1] << " skipped" << endl;
           }else{                   
              pipeline.add_dynamic_spectrum(pDynamicSpectrum);           
              std::cout << "Added dynamic spectrum for pixel " << ds_pixel[0] << "," << ds_pixel[1] << endl;
           }
        }
        std::cout << "Enabled dynamic spectrum mode.." << std::endl;
    }

    if(opts.dm_list.size() > 0){
        // enabled dedispersion
        pipeline.dedisp_engine.initialise(frequencies_list, opts.integrationTime);
    }
    for (auto& one_second_data : observation) {
        /* the std::vector memory allocator progressively allocates larger chunk of memory to
            * make space for new elements, copying existing elements to the new memory location.
            * In turn the MemoryBuffer copy contructor is called multiple times for the same object,
            * wasting time and memory. The solution is to use smart pointers to Voltages object.
            * In this way, smart pointers are the ones being copied over and over, instead of
            * the object they point to.
            */
        std::vector<std::shared_ptr<Voltages>> voltages;
        size_t ch_counter {0ull};
        high_resolution_clock::time_point read_volt_start = high_resolution_clock::now();
        #pragma omp parallel for schedule(static)
        for(size_t i = 0; i < one_second_data.size(); i++){
            auto dat_file = one_second_data[i];
            std::string& filename {dat_file.first};
            ObservationInfo obs_info {dat_file.second};
            if(opts.n_antennas > 0) obs_info.nAntennas = opts.n_antennas;
            obs_info.coarse_channel_index = i;
            unsigned int integration_steps {static_cast<unsigned int>(opts.integrationTime / obs_info.timeResolution)};
            // std::cout << "Pipeline: reading in " << filename << std::endl;
            auto volt = Voltages::from_dat_file(filename, obs_info, integration_steps);
            #pragma omp critical
            voltages.emplace_back(std::make_shared<Voltages>(std::move(volt)));
        }
        high_resolution_clock::time_point read_volt_end = high_resolution_clock::now();
        duration<double> volt_dur = duration_cast<duration<double>>(read_volt_end - read_volt_start);
        std::cout << "Reading voltages took " << volt_dur.count() << " seconds." << std::endl;
        pipeline.run(voltages);
    }
    if(pipeline.dedisp_engine.is_initialised()){
        pipeline.dedisp_engine.process_buffer();
    }
    if(dynamic_spectrum_mode_enabled) pipeline.save_dynamic_spectra();
    pipeline.finalise();
}

std::vector<float> get_frequencies(const std::vector<DatFile>& one_second, int chnls_to_avg){
    unsigned int bottom_coarse {one_second[0].second.coarseChannel};
    unsigned int n_bands = (one_second[0].second.nFrequencies / chnls_to_avg) * one_second.size();
    std::vector<float> frequencies(n_bands + 1);
    float coarse_ch_bw = one_second[0].second.coarseChannelBandwidth;
    float fine_ch_bw = one_second[0].second.frequencyResolution * chnls_to_avg;
    frequencies[0] = bottom_coarse * coarse_ch_bw - coarse_ch_bw / 2.0;
    for(int i {1}; i <= n_bands; i++){
        frequencies[i] = frequencies[i - 1] + fine_ch_bw;
    }
    for(int i {0}; i <= n_bands; i++)  frequencies[i] /= 1000.0f;
    std::cout << "Bottom frequency is " << frequencies[0] << ", top frequency is " << frequencies[frequencies.size() - 1] << std::endl;
    return frequencies;
}


void print_help(std::string exec_name, ProgramOptions& opts ){
    std::cout << "\n" << exec_name << " -t <int time> -I <input directory> [-c <cnls>] [-o <outdir>]\n"
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
    "\t-F FoV[deg]     : field of view in degrees [default 180 degree]\n" // should be opts.FOV_degrees but I am not used to cout
    "\t-w WEIGHTING    : change weighting schema N - natural, U - uniform [default N]\n" // should be opts.szWeighting    
    "\t-m MIN_UV_DISTANCE : minimum UV distance in wavelengths for a baseline to be included [default -1000]\n" // change to opts.MinUV
    "\t-a antenna_positions.txt : text file with antenna positions in a format : AntName X[m] Y[m] Z[m]\n"
    "\t-A flagged antennas list string (coma separated, e.g. 1,2,3,4) [default empty]\n"
    "\t-n IMAGE_SIZE : single value N for image size N x N pixels [default 128]\n" // change to opts.ImageSize
    "\t-Z : image phase centered at zenith, re-calculation of UVW is not required\n" 
    "\t-S : print statistics of the final sky image (mean,median,rms,rms_iqr etc)\n"
    "\t-M META_DATA_FILE : name of meta data file\n"
    "\t-v VERBOSITY/DEBUG level [default ??? ]\n" // TO-ADD CPacerImager::m_ImagerDebugLevel but I am not used to cout 
    "\t-V FILE_SAVE_LEVEL : to control number of FITS files saved [default ???]\n" // TO-ADD CPacerImager::m_SaveFilesLevel I am not used to cout 
    "\t-C frequency_channel to image [default -1 - means all]\n"
    "\t-O : oversampling factor (default: 2.0)\n"
    "'t-u : average images across frequency channels and timesteps.\n"
    "\t-P : set phase centre to RA_DEG,DEC_DEG, example -P 148.2875,7.92638889 to have B0950+08 in the phase centre\n"
    "\t-D : comma-separated list of DM trials (e.g. -D 0,0.5,1,1.5) \n"
    "\t-S : Signal to Noise (SNR) threshold level for declaring a detection.\n"
    "\t-E : polarization product to image. Options are: XX, YY, I (Stokes I). Default is I.\n"
    "\t-I : path to the directory containing the input .dat files.\n"
    "\t-X : starting second to process, expressed as number of seconds from the start of the observation. Default is 0.\n"
    "\t-Q : number of seconds of observation to process. Default is all, starting from the offset (-M).\n"
    "\t-p <postfix>: a string to optionally append to the end of output file names.\n"
    "\t-f <threshold>[,<history threshold>] enable RFI flagging by discarding all images in a "
    "FFT batch whose noise level is <threshold> times the average rms.\n"
    "\t\tOptionally, enable history thresholding by specifying a second threshold, larger than the first one. (e.g. -f5,15)\n"
    "\t-d <x,y> compute the dynamic spectrum for the (x, y) pixel.\n"
    "\t-R <n_antennas> : override the number of antennas (128)\n"
    "\t-T <cands threshold>: maximum number of candidates per batch that will be accepted (default: all candidates, no filtering)\n"
    "\t"
    << std::endl;
}



void parse_program_options(int argc, char** argv, ProgramOptions& opts){
    // default values
    opts.outputDir = std::string {"."};
    opts.nChannelsToAvg = 1;
    opts.integrationTime = -1.0; // No default!
    opts.szCalibrationSolutionsFile = "";
    opts.MinUV = -1000;
    opts.bZenithImage = false;
    opts.szWeighting = "N";
    opts.ImageSize = 180;
    opts.reorder = false;
    opts.bPrintImageStatistics = false;
    opts.FreqChannelToImage = -1; // image all channels
    opts.averageImages = false;
    opts.bChangePhaseCentre = false;
    opts.SNR = 5.0f;
    opts.pol_to_image = Polarization::I;
    opts.oversampling_factor = 2.0f;
    opts.seconds_count = -1;
    opts.seconds_offset = 0;
    opts.rfi_block_threshold = -1;
    opts.rfi_history_threshold = -1;
    opts.n_antennas = -1;
    opts.cands_per_batch = -1;
    // default debug levels :
    CPacerImager::SetFileLevel(SAVE_FILES_FINAL);
    CPacerImager::SetDebugLevel(IMAGER_WARNING_LEVEL);

    const char *options = "rt:c:o:a:M:Zi:s:F:n:v:w:V:C:A:b:uP:D:S:E:O:X:Q:I:p:f:d:R:";
    int current_opt;
    while((current_opt = getopt(argc, argv, options)) != - 1){
        switch(current_opt){
            case 'T' : {
                opts.cands_per_batch = atoi(optarg);
                break;
            }
            case 'R' : {
                opts.n_antennas = atoi(optarg);
                break;
            }
            case 'p': {
                opts.postfix = optarg;
                break;
            }
            case 'f': {
                auto items = ::tokenize_string(optarg, ',');
                if(items.size() == 1)
                    opts.rfi_block_threshold = std::stof(items[0]);
                else {
                    opts.rfi_block_threshold = std::stof(items[0]);
                    opts.rfi_history_threshold = std::stof(items[1]);
                }
                break;
            }
            case 'd': {
                auto items = ::tokenize_string(optarg, ':');
                for( auto px : items ){
                   auto px_coords = ::tokenize_string(px, ','); 
                   if(px_coords.size() != 2) throw std::invalid_argument {"Invalid pixel specification for the -d option."};
                   
                   int x = atoi(px_coords[0].c_str());
                   int y = atoi(px_coords[1].c_str());
                   opts.ds_pixels.push_back(array<int,2>{x,y});
                }
                break;
            }
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

            case 'D': {
               if( optarg && strlen(optarg)){
                  opts.dm_list_string = optarg;
               }
               break;
            }
            
            case 'b' : {
               opts.coarseChannelIndex = atoi(optarg);
               break;
            }
            case 'S' : {
                opts.SNR = atof(optarg);
                break;
            }

            case 'C': {
               if( optarg && strlen(optarg) ){
                  opts.FreqChannelToImage = atol( optarg );
               }
               break;
            }
            case 'O': {
                opts.oversampling_factor = atof(optarg);
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
            case 'F': {
                opts.FOV_degrees = atof(optarg);
                break;
            }
            case 'n': {
                opts.ImageSize = atoi(optarg);
                break;
            }
            case 'X': {
                opts.seconds_offset = atoi(optarg);
                break;
            }
            case 'Q': {
                opts.seconds_count = atoi(optarg);
                break;
            }
            case 'I': {
                opts.input_directory = std::string {optarg};
                break;
            }
            case 'P' : {
               if(optarg && strlen(optarg)){
                    std::string radec_str {optarg};
                   parse_radec_string(radec_str, opts.ra_dec);
                   opts.bChangePhaseCentre = true;
                   printf("DEBUG : changing phase centre to (RA,DEC) = (%.8f,%.8f) [deg]\n",opts.ra_dec[0], opts.ra_dec[1]);
               }
               break;
            }
            case 'E': {
                if(optarg == std::string {"XX"}){
                    opts.pol_to_image = Polarization::XX;
                }else if(optarg == std::string {"YY"}){
                    opts.pol_to_image = Polarization::YY;
                }else if(optarg == std::string {"I"}){
                    opts.pol_to_image = Polarization::I;
                }else{
                    throw std::invalid_argument("Value for polarization to image not recognised. Valid options are 'XX', 'YY', and 'I'.");
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

            case 'u': {
                opts.averageImages = true;
                break;
            }

            default : {
                std::stringstream ss;
                ss << "Unrecognised option: '" << static_cast<char>(optopt) << "'.";
                throw std::invalid_argument(ss.str());
            }
        }
    }
    
    // options validation
    if(opts.integrationTime <= 0) throw std::invalid_argument("You must specify a value for the integration time.");
    if(opts.input_directory == "") throw std::invalid_argument("No input directory specified.");
    if(!blink::dir_exists(opts.input_directory)) throw std::invalid_argument("Invalid input directory specified.");

    if(strlen(opts.gFlaggedAntennasListString.c_str()) > 0 ){
        auto items = ::tokenize_string(opts.gFlaggedAntennasListString);
       for(int i=0;i<items.size();i++){
          opts.gFlaggedAntennasList.push_back(atol( items[i].c_str()));
       }
    }

    if(opts.dm_list_string.length() > 0 ){
       if(opts.dm_list_string.find(":") != std::string::npos){
            // dm trials specified with min:max:step
            auto items = ::tokenize_string(opts.dm_list_string, ':');
            if(items.size() != 3) throw std::invalid_argument("Invalid DM range specified.");
            float dm_start = atof(items[0].c_str());
            float dm_stop = atof(items[1].c_str());
            float dm_delta = atof(items[2].c_str());
            while(dm_start <= dm_stop){
                std::cout << "DM trial " << dm_start << std::endl;
                opts.dm_list.push_back(dm_start);
                dm_start += dm_delta;
            }
       }else{
            // dm trials specified explicitly with a list
            auto items = ::tokenize_string(opts.dm_list_string);
            for(int i=0;i<items.size();i++){
                opts.dm_list.push_back(atof(items[i].c_str()));
            }
       }
    }
}


void print_program_options(const ProgramOptions& opts){
    std::cout << "\nRunning the correlator program with the following options:\n"
    "\t Integration time interval: " << opts.integrationTime << "s\n"
    "\t Freq. channel: " << opts.FreqChannelToImage << "\n"
    "\t Number of channels to average: " << opts.nChannelsToAvg << "\n"
    "\t Output directory: " << opts.outputDir << "\n"
    "\t Calibration file: " << opts.szCalibrationSolutionsFile << "\n"    
    "\t Zenith image: " << opts.bZenithImage << std::endl;
}



double ra_to_degrees(const std::string& ra_str) {
    // RA format: HH:MM:SS(.ss), where HH=[00–23], MM=[00–59], SS=[00–59.99]
    std::regex ra_pattern(R"(^(2[0-3]|1\d|0?\d):([0-5]?\d):([0-5]?\d(?:\.\d{1,2})?)$)");
    std::smatch match;
    if (!std::regex_match(ra_str, match, ra_pattern)) {
        throw std::invalid_argument("Invalid RA format: " + ra_str);
    }

    int hours   = std::stoi(match[1]);
    int minutes = std::stoi(match[2]);
    double seconds = std::stod(match[3]);

    return (hours + minutes / 60.0 + seconds / 3600.0) * 15.0;
}

double dec_to_degrees(const std::string& dec_str) {
    // DEC must match ±DD:MM:SS(.ss), DD in [0,90], MM in [0,59], SS in [0,59.99]
    std::regex dec_pattern(R"(^([+-])(\d{1,2}):(\d{1,2}):(\d{1,2}(?:\.\d{1,2})?)$)");
    std::smatch match;
    if (!std::regex_match(dec_str, match, dec_pattern)) {
        throw std::invalid_argument("Invalid DEC format: " + dec_str);
    }

    char sign = match[1].str()[0];
    int degrees = std::stoi(match[2]);
    int minutes = std::stoi(match[3]);
    double seconds = std::stod(match[4]);

    double value = degrees + minutes / 60.0 + seconds / 3600.0;
    return (sign == '-') ? -value : value;
}

void parse_radec_string(const std::string& radec_string, double ra_dec[]){
    auto ra_dec_component_str = ::tokenize_string(radec_string, ',');
    if(ra_dec_component_str.size() != 2)
        throw std::invalid_argument {"Invalid RA DEC specification."};
    if(radec_string.find(":") != std::string::npos){
        // hours format
        ra_dec[0] = ra_to_degrees(ra_dec_component_str[0]);
        ra_dec[1] = dec_to_degrees(ra_dec_component_str[1]);
    }else{
        // already in degree format
        ra_dec[0] = std::stod(ra_dec_component_str[0]);
        ra_dec[1] = std::stod(ra_dec_component_str[1]);
    }
}
