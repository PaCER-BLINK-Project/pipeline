// program averages few FITS images of the same sizes 

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <unistd.h>

#include <vector>
#include <iostream>
#include <fstream> 
#include <sstream>
using namespace std;

string infile="total_power.txt";
string outfile="exclude_ranges.txt";
double gThreshold=5.00;

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
      cout << "Open OK !!!";
   else 
      cerr << "ERROR : could not open file " << filename;

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

      printf("LINE : |%s|\n",line.c_str());
      printf("     : |%ld %.4f %.4f %.4f %.4f|\n",(total_power_vec.end()-1)->timeindex,(total_power_vec.end()-1)->total_power,(total_power_vec.end()-1)->median_of_median,(total_power_vec.end()-1)->median_of_rmsiqr,(total_power_vec.end()-1)->median);
//      printf("     : |%ld %.5f %.5f %.5f %.5f|\n",tmp.timeindex,tmp.total_power,tmp.median_of_median,tmp.median_of_rmsiqr,tmp.median);
//      printf("     : |%s %s %s %s %s|\n",items[0].c_str(),items[1].c_str(),items[2].c_str(),items[3].c_str(),items[4].c_str());
   }
   
   ifile.close();
      
   return total_power_vec.size();
}

void find_exclude_ranges( vector<cTotalPower>& total_power_vec, const char* outfile ) 
{
   int n_points = total_power_vec.size();
   FILE* outf = fopen( outfile ,"w");
   printf("DEBUG : n_points = %d\n",n_points);
   printf("Exclusion ranges:\n");
   // generate exclusion ranges :
   int border = 10;
   int start_range = -1;
   int last_range_start = -1, last_range_end = -1;
   int count_exclude_range=0;
   for(int i=50;i<n_points;i++){  // skips first 50 (number of points in the running median) - so before running median / iqr stabilise 
      if( total_power_vec[i].total_power > total_power_vec[i].up ){
         if( start_range < 0 ){
            if( last_range_end>=0 && (i-last_range_end)<10 ){
               // if small break between ranges - extend the current range :
               start_range = last_range_start;
            }else{
               if( last_range_start>=0 && last_range_end>=0 ){
                  // save current range :
                  last_range_start -= border; // a bit agressive border
                  last_range_end -= border;
                  printf("\tExclusion range = %d - %d\n",last_range_start,last_range_end);
                  fprintf(outf,"%d-%d\n",last_range_start,last_range_end);
                  last_range_start = -1;
                  last_range_end = -1;
                  count_exclude_range++;
               }
               start_range = i; 
            }
         }
      }else{
         if( start_range >= 0 && (i-start_range)>10 ){
            int t_start = total_power_vec[start_range].timeindex - border;
            int t_end = total_power_vec[i-1].timeindex + border;
//            printf("\tExclusion range = %d - %d\n",t_start,t_end);
//            fprintf(outf,"%d-%d\n",t_start,t_end);
            last_range_start = t_start;
            last_range_end = t_end;
            start_range = -1;
         }
      }

      if( last_range_start>=0 && last_range_end>=0 && (i-last_range_end)>=10 ){
          // save current range :
          last_range_start -= border;
          last_range_end -= border;
          printf("\tExclusion range = %d - %d\n",last_range_start,last_range_end);
          fprintf(outf,"%d-%d\n",last_range_start,last_range_end);
          last_range_start = -1;
          last_range_end = -1;
          count_exclude_range++;
      }

//      printf("%d : %.8f vs. %.8f\n",i,total_power[i],up[i]);
   }
   fclose(outf);

   printf("Number of exclusion ranges = %d\n",count_exclude_range);

}

void usage()
{
   printf("exclude_ranges_totalpower TOTAL_POWER.txt\n\n\n");
   printf("TOTAL_POWER.txt file should be the output from test_totalpower program with 5 columns : # TIMEINDEX TOTAL_POWER MedianOfMedians MedianOfRMSIQR MEDIAN\n");
   printf("-o OUTPUT_FILE_NAME [default %s]\n",outfile.c_str());
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
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }
}

void print_parameters()
{
    printf("############################################################################################\n");
    printf("PARAMETERS :\n");
    printf("############################################################################################\n");
    printf("Input file    = %s\n",infile.c_str());
    printf("Output file   = %s\n",outfile.c_str());
    printf("############################################################################################\n");
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
  
  find_exclude_ranges( total_power_vec, outfile.c_str() );
}

