#ifndef __PIPELINE__
#define __PIPELINE__

#include <vector>
#include <string>
#include <astroio.hpp>
#include <utils.hpp>
#include <memory>
#ifdef IMAGER_HIP
#include <hip/pacer_imager_hip.h>
#else
#include <pacer_imager.h>
#endif

namespace blink {
    
    enum DataType {MWA=1, EDA2=2};

    struct ProgramOptions {
        std::vector<std::string> input_files;
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
        // all-sky EDA2-like image
        DataType inputDataType;
        int FreqChannelToImage;
        
        // flagging antennas:
        string gFlaggedAntennasListString;
        vector<int> gFlaggedAntennasList;
    };


    class Pipeline {

        #ifdef IMAGER_HIP
           CPacerImagerHip imager;   
        #else
           CPacerImager imager;
        #endif  

        std::string output_dir;
        bool calibrate {false};        
        bool reorder {false};
        // Correlation-related options
        unsigned int channels_to_avg {1};
        double integration_time {0.01};

        
        // might be better to give it as input?
        std::string calibration_solutions_file;

        int imageSize {512};
        // might be better to give it as input
        std::string MetaDataFile;
        

        std::string ImagerOutFilePostfix;
        double MinUV;
        double frequencyMHz;
        double FOV_degrees;
        bool bPrintImageStatistics;
        std::string szWeighting;
        bool bZenithImage; // all-sky EDA2-like image
        
        DataType inputDataType;
        
        // flagged antennas :
        std::string szFlaggedAntennasListString;
        std::vector<int> szFlaggedAntennasList;
        

        public:        

        Pipeline(unsigned int nChannelsToAvg, double integrationTime, bool reorder, bool calibrate, std::string solutions_file,
                  int imageSize, std::string metadataFile, std::string szAntennaPositionsFile, double minUV, 
                  bool printImageStats, std::string szWeighting, std::string outputDir, bool bZenithImage,
                  double FOV_degrees, blink::DataType inputType,
                  vector<int>& flagged_antennas, std::string& output_dir
                );
        
        void run(const Voltages& input, int freq_channel = -1);
        void run(const std::vector<std::shared_ptr<Voltages>>& inputs, int freq_channel = -1);

    };
}


#endif
