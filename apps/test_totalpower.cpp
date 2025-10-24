// program averages few FITS images of the same sizes 

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include <bg_globals.h>
#include <bg_fits.h>
#include "../src/totalpower.h"
#include "../src/dynamic_spectrum.hpp"

#include <vector>
using namespace std;

string list="fits_list";

void usage()
{
   printf("test_totalpower fits_list\n\n\n");
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "h";
   int opt;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'h':
            usage();
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
    printf("List file    = %s\n",list.c_str());
    printf("############################################################################################\n");
}

int main(int argc,char* argv[])
{
  string list="fits_list";
  if( argc >= 2 ){
     list = argv[1];
  }
  parse_cmdline( argc , argv );
  print_parameters();
 
  
  vector<string> fits_list;
  if( bg_read_list(list.c_str(),fits_list) <= 0 ){
     printf("ERROR : could not read list file %s\n",list.c_str());
     exit(-1);
  }else{
     for(int i=0;i<fits_list.size();i++){
        printf("%i %s\n",i,fits_list[i].c_str());
     }
  }

  CBgFits first_fits;
  if( first_fits.ReadFits( fits_list[0].c_str() , 0, 1, 1 ) ){
     printf("ERROR : could not read fits file %s on the list\n",fits_list[0].c_str());
     exit(-1); 
  }
  
  HeaderRecord* pKey = first_fits.GetKeyword( "CDELT2" );
  double cdelt1 = 20.00;
  if( pKey ){
     cdelt1 = atof( pKey->Value.c_str() );
  }

  DynamicSpectrum dynaspec( first_fits.GetXSize() , first_fits.GetYSize(), 1, 0, 0  );

  TotalPower total_power_history(cdelt1);

  int ntimes=100;
  long int start_time = 0;  
  for(int i=0;i<fits_list.size();i++){
     CBgFits fits;
     if( fits.ReadFits( fits_list[i].c_str() , 0, 1, 1 ) ){
        printf("ERROR : could not read fits file %s on the list\n",fits_list[i].c_str());
        exit(-1); 
     }else{
        printf("OK : fits file %s read ok\n",fits_list[i].c_str());
     }

// FULL FITS VERSION :     
     total_power_history.filename = fits_list[i].c_str();
     total_power_history.filename.replace(total_power_history.filename.find(".fits"),5,".total_power");
     total_power_history.ResetFile();
     printf("DEBUG : saving total power from %s to %s\n",fits_list[i].c_str(),total_power_history.filename.c_str());
     
//     total_power_history.calc( fits, true, true, 0, -1, start_time );
     int blocks = int( fits.GetXSize()/ntimes );
     int offset = 0; 
     for( int b=0;b<blocks;b++){
        total_power_history.calc( fits, true, true, offset, ntimes, start_time );
        offset += ntimes;
     }
//     start_time += fits.GetXSize();          
  }
}

