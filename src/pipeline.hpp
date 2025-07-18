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
        
        
        // flagged antennas :
        std::string szFlaggedAntennasListString;
        std::vector<int> szFlaggedAntennasList;
        
        // dedispersion parameters
        int _current_top_frequency_band_idx {769};

        int sweep_size;
        int batch_size;
        int buffer_size;
        int table_size;

        int window_start_idx {0};
        // start position of the current batch within the window
        int window_offset {0};

        std::vector<float> frequencies;
        std::vector<float> dm_list;
        std::vector<int> delay_table;
        MemoryBuffer<float> dm_starttime;
        

        public:        

        Pipeline(unsigned int nChannelsToAvg, double integrationTime, bool reorder, bool calibrate, std::string solutions_file,
                  int imageSize, std::string metadataFile, std::string szAntennaPositionsFile, double minUV, 
                  bool printImageStats, std::string szWeighting, std::string outputDir, bool bZenithImage,
                  double FOV_degrees, bool averageImages,
                  vector<int>& flagged_antennas, std::vector<float>& dm_list, std::string& output_dir
                );
        
        void run(const Voltages& input);
        void run(const std::vector<std::shared_ptr<Voltages>>& inputs);

        void set_frequencies(const std::vector<float>& frequencies);

    };
}


#endif
