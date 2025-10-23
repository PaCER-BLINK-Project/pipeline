#ifndef TOTAL_POWER_H_
#define TOTAL_POWER_H_

#include <vector>
#include <deque>
#include <string>

// single record for a 1-second block of VCS data (calculated based on the full time resolution data)
struct MedianIQRRecord
{
   double median;
   double rmsiqr;
};

// history of medians, RMSIQRs per 1-second block of VCS data 
class MedianIQRHistory : public std::deque<MedianIQRRecord>
{
   std::vector<double> sorting_table;
public:
   MedianIQRHistory( int maxElements );
   ~MedianIQRHistory();
   void add( double median , double rmsiqr );
   
   double m_MedianOfMedians; // median of median values in history, updated on each add
   double m_MedianOfRmsIQRs; // median of rmsiqr values in history, updated on each add

   int m_MaxElements;
};

// total power record in full time resolution
struct TotalPowerRec
{
public:     
   TotalPowerRec();

   long int time_index;
   double   total_power;
};

class DynamicSpectrum;
class CBgFits;

class TotalPower : std::vector<TotalPowerRec> {

public:
   std::string filename {"total_power.txt"} ; // output file name 
   std::vector<TotalPowerRec> m_PreviousBuffer; // total power for the previous 1-second block 
  
   TotalPower( double time_resolution_ms, int history_size=100, double start_freq_mhz=130, double end_freq_mhz=170 ); // history of 100 previous seconds - as in Fredda version
   ~TotalPower();
   void ResetFile();   

   // find maximum total power in some radius    
   double get_max_total_power(int index, int radius);
   double get_max_total_power_before(int index, int radius);
   double get_max_total_power_after(int index, int radius);

   double dispersion_delay_ms(double dm);
 
   double  m_StartFreqMHz;
   double  m_EndFreqMHz;
   double  m_TimeResolutionMS;
   std::vector<double> sorted_total_power;
   double m_Median;
   double m_RMSIQR;
   
   double GetMedianOfMedians();
   double GetMedianOfRMSIQRs();

   MedianIQRHistory m_MedianIQRHistory;

// Public interface : 
   void calc( float *dynaspec, size_t n_timesteps, size_t n_channels, bool do_dump=false, bool use_rms=true, int offset=0, int ntimes=-1, long int start_time=0 );
   void calc( DynamicSpectrum& dynaspec, bool do_dump=false, bool use_rms=true, int offset=0, int ntimes=-1, long int start_time=0 );
   void calc( CBgFits& dynaspec, bool do_dump=false, bool use_rms=true, int offset=0, int ntimes=-1, long int start_time=0 );
   
   bool is_total_power_ok(DynamicSpectrum& dynaspec, int t, double total_power_threshold=5 );
};

#endif

