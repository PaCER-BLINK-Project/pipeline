#include <stdio.h>
#include <cassert>
#include "totalpower.h"
#include "dynamic_spectrum.hpp"
#include <bg_fits.h>

#include <algorithm>

MedianIQRHistory::MedianIQRHistory( int maxElements )
: m_MaxElements(maxElements), m_MedianOfMedians(-1), m_MedianOfRmsIQRs(-1)
{
}

MedianIQRHistory::~MedianIQRHistory()
{
}

void MedianIQRHistory::add( double median , double rmsiqr )
{
   MedianIQRRecord tmp;
   tmp.median = median;
   tmp.rmsiqr = rmsiqr;

   printf("DEBUG_MedianIQRHistory: ");
   for(auto it = cbegin();it!=cend();it++){
     printf("%.3f ",it->median);
   }
   printf("\n");
   if( size() >= m_MaxElements ){
      pop_front(); // remove element in front
   }
   push_back( tmp );

   assert( size() <= m_MaxElements );

   // Sort median and rmsiqr values to get median of medians and median of RMSQRSs
   // TODO: ? probably there could be some more efficient and clever way of "insertion" like into a binary tree.
   //       or implement fit or Savitzky-Goley or perfect smoothing (see Adrian S.)
   // once it works, I can implement improvements :
   //    - keep array of values in double* array and list of indexes in queue
   //    - add new element : remove index from the front, push to back of queue of
   //                        indexes , replaces in array with new element,
   //                        move the element like in bubble sort to find its place , update
   //                        indexes in the queue would be the most difficult part in this way
   // sort Median values :
   sorting_table.clear();
   for(auto it = cbegin();it!=cend();it++){
      sorting_table.push_back(it->median);
   }
   ::sort(sorting_table.begin(),sorting_table.end());
   
   printf("DEBUG_MEDIAN_OF_MEDIANS:\n");
   for(int k=0;k<sorting_table.size();k++){
      printf("%.3f ",sorting_table[k]);
   }
   m_MedianOfMedians = sorting_table[int(sorting_table.size()/2)];
   printf("\n -> median_of_medians = %.3f\n",m_MedianOfMedians);
   
   sorting_table.clear();
   for(auto it = cbegin();it!=cend();it++){
      sorting_table.push_back(it->rmsiqr);
   }
   ::sort(sorting_table.begin(),sorting_table.end());

   
   m_MedianOfRmsIQRs = sorting_table[int(sorting_table.size()/2)];     
}

TotalPowerRec::TotalPowerRec()
: time_index(-1), total_power(0.00)
{
}

TotalPower::TotalPower( double time_resolution_ms, int history_size, double start_freq_mhz, double end_freq_mhz )
: m_Median(0.00), m_RMSIQR(0.00), m_StartFreqMHz(start_freq_mhz), m_EndFreqMHz(end_freq_mhz), m_TimeResolutionMS(time_resolution_ms), m_MedianIQRHistory(history_size)
{
   printf("DEBUG : created TotalPower object in %.4f - %.4f MHz range, time resolution = %.4f [ms], will keep %d samples history\n",start_freq_mhz,end_freq_mhz,m_TimeResolutionMS,history_size);
}

TotalPower::~TotalPower()
{
}


double TotalPower::GetMedianOfMedians()
{
   return m_MedianIQRHistory.m_MedianOfMedians;
}

double TotalPower::GetMedianOfRMSIQRs()
{
   return m_MedianIQRHistory.m_MedianOfRmsIQRs;
}

void TotalPower::ResetFile()
{
   // just clean up file 
   FILE* out_f = fopen(filename.c_str(),"w");
   // rec.time_index,rec.total_power,GetMedianOfMedians(),GetMedianOfRMSIQRs(),m_Median);
   fprintf(out_f,"# TIMEINDEX TOTAL_POWER MedianOfMedians MedianOfRMSIQR MEDIAN\n");
   fclose(out_f);
}


void TotalPower::calc( float *dynaspec, size_t n_timesteps, size_t n_channels, bool do_dump, bool use_rms, int offset, int ntimes, long int start_time )
{
   if( ntimes <= 0 ){
      ntimes = n_timesteps;
   }
   int nchan  = n_channels;
 
   // copy current vector to m_PreviousBuffer to still keep it for checking
   // early samples in the new buffer 
   if( size() > 0 ){
      m_PreviousBuffer.clear();  
      copy(begin(),end(),back_inserter(m_PreviousBuffer));
      
      printf("DEBUG : added %d records to previous total power buffer\n",int(size()));
   }

   clear();

   TotalPowerRec tmp;
   int time_size = n_timesteps;
   float* data = dynaspec;
   printf("DEBUG : TotalPower::calc, ntimes = %d, nchan = %d, first values %.2f, %.2f, %.2f\n",ntimes,nchan,data[0],data[1],data[2]);
  
   for(long int t=offset;t<(offset+ntimes);t++){
      double sum = 0.00;
      double sum2 = 0.0;
      int count = 0; // in case some excluded in the future and nchan!=count 
      for(int ch=0;ch<nchan;ch++){
        double value = data[time_size*ch+t];
        sum += value;
        sum2 += value*value;
        count++;
        if(t==0 && false){printf("\t%.2f\n",value);}
      }
      if(t==0 && false){printf("\n");}

      tmp.time_index = start_time + t;;
      double rms = sqrt( sum2/count - (sum/count)*(sum/count) );
      if( use_rms )
         tmp.total_power = rms;
      else 
         tmp.total_power = sum/nchan; // just to keep consistency and compare to dumpfilfile_float (main_float.cpp in mwafrb/src)


      if(t==0){printf("First total power = %.8f\n",tmp.total_power);}
      push_back(tmp);
   }   

   sorted_total_power.clear();
   for( auto power : (*this) ){
      sorted_total_power.push_back(power.total_power);
   }
   ::sort( sorted_total_power.begin(), sorted_total_power.end());

   m_Median = sorted_total_power[int(sorted_total_power.size()/2)];
   int q75= int(sorted_total_power.size()*0.75);
   int q25= int(sorted_total_power.size()*0.25);
   m_RMSIQR = ( sorted_total_power[q75] - sorted_total_power[q25] ) / 1.35;

   // add to history used for getting filtered / smoothed value (stable and excluding spikes)
   m_MedianIQRHistory.add( m_Median , m_RMSIQR );

   printf("Double check of first total power = %.8f, median of MEDIANS = %.4f (vs. local median = %.4f), Median of RMSIQRs = %.4f (vs. local rmsiqr = %.4f )\n",(*this)[0].total_power,GetMedianOfMedians(),m_Median,GetMedianOfRMSIQRs(),m_RMSIQR);
   if( do_dump ){
      FILE* out_f = fopen(filename.c_str(),"a+");
      for(int i=0;i<size();i++){
         TotalPowerRec& rec = (*this)[i];
         fprintf(out_f,"%ld %.4f %.8f %.8f %.8f\n",rec.time_index,rec.total_power,GetMedianOfMedians(),GetMedianOfRMSIQRs(),m_Median);
      }
      fclose(out_f);
   }
}

void TotalPower::calc( DynamicSpectrum& dynaspec, bool do_dump, bool use_rms, int offset, int ntimes, long int start_time )
{
   calc( dynaspec.data(), dynaspec.get_ntimesteps(), dynaspec.get_nchannels(), do_dump, use_rms, offset, ntimes, start_time );
}

void TotalPower::calc( CBgFits& dynaspec, bool do_dump, bool use_rms, int offset, int ntimes, long int start_time )
{
    calc( dynaspec.get_data(), dynaspec.GetXSize(), dynaspec.GetYSize(), do_dump, use_rms, offset, ntimes, start_time );
}

bool TotalPower::is_total_power_ok(DynamicSpectrum& dynaspec, int t, double total_power_threshold )
{
   if( t>=0 && t<=size() ){
      double total_power = (*this)[t].total_power;
      double median = GetMedianOfMedians();
      double local_median = m_Median;
      double rmsiqr = GetMedianOfRMSIQRs();
      
      double threshold = total_power_threshold*rmsiqr + median;
      
      if( total_power > threshold ){
         return false;
      }
   }
   
   return true;
}

double TotalPower::get_max_total_power(int index, int radius)
{
   double max_power = -1;
   int start = index - radius;

   for(int i=start;i<(index+radius);i++){
      if( i>=0 && i<size() ){
         TotalPowerRec& rec = (*this)[i];
         if( rec.total_power > max_power ){
            max_power = rec.total_power;
         } 
      }
   }

   return max_power;
}

double TotalPower::get_max_total_power_before(int index, int radius)
{
   double max_power = -1;
   int start = index - radius;

   for(int i=start;i<(index);i++){
      if( i>=0 && i<size() ){
         TotalPowerRec& rec = (*this)[i];
         if( rec.total_power > max_power ){
            max_power = rec.total_power;
         } 
      }else{
         // if i<0 -> try to check in the total power from the previous buffer copied to m_PreviousBuffer
         if( i < 0 && m_PreviousBuffer.size()>0 ){
            int prev_buffer_i = m_PreviousBuffer.size() + i;
            if( prev_buffer_i>=0 && prev_buffer_i<m_PreviousBuffer.size() ){
               TotalPowerRec& rec = m_PreviousBuffer[prev_buffer_i];
               if( rec.total_power > max_power ){
                  max_power = rec.total_power;
               }
            }
         }
      }
   }

   return max_power;

}

double TotalPower::get_max_total_power_after(int index, int radius)
{
   double max_power = -1;
   int start = index;

   for(int i=start;i<(index+radius);i++){
      if( i>=0 && i<size() ){
         TotalPowerRec& rec = (*this)[i];
         if( rec.total_power > max_power ){
            max_power = rec.total_power;
         } 
      }
   }

   return max_power;  
}


double TotalPower::dispersion_delay_ms(double dm)
{
   double freq2_ghz = m_EndFreqMHz / 1000.00; 
   double freq1_ghz = m_StartFreqMHz / 1000.00;

   double delta_t_ms = 4.15 * dm * ( 1.00/(freq1_ghz*freq1_ghz) - 1.00/(freq2_ghz*freq2_ghz) );

   printf("DEBUG : delta_t_ms = %.6f [ms] for %.4f - %.4f MHz range, and dm = %.3f [ps/cm^3]\n",delta_t_ms,freq1_ghz,freq2_ghz,dm);

   return delta_t_ms;
}
