// program averages few FITS images of the same sizes 

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <unistd.h>
#include <iomanip>

#include <vector>
#include <iostream>
#include <fstream> 
#include <sstream>
using namespace std;

string infile="total_power.txt";
string outfile="exclude_ranges.txt";
double gThreshold=5.00;
bool   gDebug=false;
int    gRunningMedianSize=50;
int    gBorder=10;

struct cTotalPower
{
   long int timeindex;
   double   total_power;
   double   median_of_median;
   double   median_of_rmsiqr;
   double   median;  
   
   double up;
   double down;
};

int read_file( const char* filename, vector<cTotalPower>& total_power_vec )
{
   total_power_vec.clear();
   
   ifstream ifile{filename};
   if(ifile)
      std::cout << "Open OK !!!" << std::endl;
   else 
      std::cerr << "ERROR : could not open file " << filename << std::endl;

   string line;      
   while(getline(ifile,line)){ // returns false when file ends, removes \n from the text (no newline)
//      cout << "DEBUG: " << line << endl;
      
      std::vector<std::string> items;
      std::string item;

      // Create an input string stream from the data string
      std::istringstream iss(line);

      // Extract items using getline with a comma delimiter
      while (std::getline(iss, item, ' ')) {
          items.push_back(item);
      }       
      if( items.size() < 5 ){      
         cout << "WARNING : line has " << items.size() << " items -> less than 5 items -> ignored" << endl;
         continue;
      }
      if( items[0] == "#" ){
//         cout << "COMMENT LINE SKIPPED" << endl;
         continue;
      }
 
      long int timeindex = atol(items[0].c_str());
      double total_power = atof(items[1].c_str());
      double median_of_medians = atof(items[2].c_str());
      double median_of_rmsiqr  = atof(items[3].c_str());
      double median = atof(items[4].c_str());
      double up = median_of_medians + gThreshold*median_of_rmsiqr;
      double down = median_of_medians - gThreshold*median_of_rmsiqr;
      total_power_vec.push_back( cTotalPower{ timeindex, total_power, median_of_medians, median_of_rmsiqr, median, up, down } );

      if( gDebug ){
         std::cout << "LINE : |" << line << "|" << std::endl;
         std::cout << std::setprecision(4);
         std::cout << "       |" << (total_power_vec.end()-1)->timeindex << " " << (total_power_vec.end()-1)->total_power << " " << (total_power_vec.end()-1)->median_of_median << " " << (total_power_vec.end()-1)->median_of_rmsiqr << " " << (total_power_vec.end()-1)->median << std::endl;
      }
   }
   
   ifile.close();
      
   return total_power_vec.size();
}

void find_exclude_ranges( vector<cTotalPower>& total_power_vec, const char* outfile, int running_median_size, int border ) 
{
   int n_points = total_power_vec.size();
   ofstream ofile(outfile);
   std::cout << "DEBUG : n_points = " << n_points << std::endl;
   std::cout << "Exclusion ranges:" << std::endl;
   // generate exclusion ranges :
   int start_range = -1;
   int last_range_start = -1, last_range_end = -1;
   int count_exclude_range=0;
   for(int i=running_median_size;i<n_points;i++){  // skips first 50 (number of points in the running median) - so before running median / iqr stabilise 
      if( total_power_vec[i].total_power > total_power_vec[i].up ){
         if( start_range < 0 ){
            if( last_range_end>=0 && (i-last_range_end)<border ){
               // if small break between ranges - extend the current range :
               start_range = last_range_start;
            }else{
               if( last_range_start>=0 && last_range_end>=0 ){
                  // save current range :
                  last_range_start -= border; // a bit agressive border
                  last_range_end -= border;
                  std::cout << "\tExclusion range = " << last_range_start << " - " << last_range_end << std::endl;
                  ofile << last_range_start << "-" << last_range_end << std::endl;
                  last_range_start = -1;
                  last_range_end = -1;
                  count_exclude_range++;
               }
               start_range = i; 
            }
         }
      }else{
         if( start_range >= 0 && (i-start_range)>border ){
            int t_start = total_power_vec[start_range].timeindex - border;
            int t_end = total_power_vec[i-1].timeindex + border;
            last_range_start = t_start;
            last_range_end = t_end;
            start_range = -1;
         }
      }

      if( last_range_start>=0 && last_range_end>=0 && (i-last_range_end)>=border ){
          // save current range :
          last_range_start -= border;
          last_range_end -= border;
          std::cout << "\tExclusion range = " << last_range_start << " - " << last_range_end << std::endl;
          ofile << last_range_start << "-" << last_range_end << std::endl;
          last_range_start = -1;
          last_range_end = -1;
          count_exclude_range++;
      }

   }
   ofile.close();

   std::cout << "Number of exclusion ranges = " << count_exclude_range << std::endl;

}

void usage()
{
   std::cout << "exclude_ranges_totalpower TOTAL_POWER.txt" << std::endl << std::endl << std::endl;   
   std::cout << "TOTAL_POWER.txt file should be the output from test_totalpower program with 5 columns : # TIMEINDEX TOTAL_POWER MedianOfMedians MedianOfRMSIQR MEDIAN" << std::endl;
   std::cout << "-o OUTPUT_FILE_NAME [default " << outfile << "]" << std::endl;
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "ho:";
   int opt;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'h':
            usage();
            break;

         case 'o':
            if( optarg ){
               outfile = optarg;
            }
            break;


         default:   
            std::cerr << "Unknown option " << opt << std::endl;
            usage();
      }
   }
}

void print_parameters()
{
    std::cout << "############################################################################################" << std::endl;
    std::cout << "PARAMETERS :" << std::endl;
    std::cout << "############################################################################################" << std::endl;
    std::cout << "Input file    = %s" << infile << std::endl;
    std::cout << "Output file   = %s" << outfile  << std::endl;
    std::cout << "############################################################################################" << std::endl;
}

int main(int argc,char* argv[])
{
  if( argc >= 2 ){
     infile = argv[1];
  }
  parse_cmdline( argc , argv );
  print_parameters();
 
  vector<cTotalPower> total_power_vec;  
  read_file( infile.c_str() , total_power_vec );
  
  find_exclude_ranges( total_power_vec, outfile.c_str(), gRunningMedianSize, gBorder );
}

