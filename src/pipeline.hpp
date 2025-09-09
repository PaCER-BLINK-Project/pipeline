#ifndef __PIPELINE__
#define __PIPELINE__

#include <vector>
#include <string>
#include <astroio.hpp>
#include <utils.hpp>
#include <memory>
#include <calibration.hpp>
#include <gpu/pacer_imager_hip.h>
#include <dedispersion.hpp>

using namespace blink::dedispersion;

namespace blink {
    
    enum DataType {MWA=1, EDA2=2};

    class Pipeline {

        std::vector<CPacerImagerHip*> imager;
        std::string output_dir, postfix;
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
        std::vector<MemoryBuffer<int>> mapping;
        std::vector<CalibrationSolutions> cal_sol;


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
        // RFI flagging threshold. A negative value disables
        // flagging.
        float rfi_flagging {-1.0};
        int num_gpus;
        
        
        public:        
        Dedispersion dedisp_engine;

        Pipeline(unsigned int nChannelsToAvg, double integrationTime, bool reorder, bool calibrate, std::string solutions_file,
                  int imageSize, std::string metadataFile, float oversampling_factor, std::string szAntennaPositionsFile,
                  double minUV, bool printImageStats, std::string szWeighting, std::string outputDir, bool bZenithImage,
                  double FOV_degrees, bool averageImages, Polarization pol_to_image,
                  vector<int>& flagged_antennas,bool change_phase_centre, double ra_deg, double dec_deg,
                  Dedispersion& dedisp_engine, float rfi_flagging, std::string& output_dir, std::string& postfix
                );
        
        void run(const Voltages& input, int gpu_id);
        void run(const std::vector<std::shared_ptr<Voltages>>& inputs);

        ~Pipeline() {
            for(CPacerImagerHip* p : imager) delete p;
        };

    };
}


#endif
