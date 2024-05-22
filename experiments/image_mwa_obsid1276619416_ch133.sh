#!/bin/bash

#SBATCH --account=director2183
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120gb
#SBATCH --output=./image_mwa_obsid1276619416_ch133.o%j
#SBATCH --error=./image_mwa_obsid1276619416_ch133.e%j
#SBATCH --export=NONE

# /home/msok/Desktop/PAWSEY/PaCER/logbook/20230908_testing_BLINK_PIPELINE_cable_and_geom_correction_applycal.odt
# see /media/msok/80f59a5b-8bba-4392-8b31-a4f02fbee2f1/mwa/vcs/1276619416/1276619416_new/cotter/cotter_only_geomandcablecorr_applycal/BLINK_PIPELINE

module_path=`which module`

if [[ $PAWSEY_CLUSTER = "setonix" ]]; then
   module reset 
   module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0 /software/setonix/2023.08/modules/zen3/gcc/12.2.0/libraries
   module load msfitslib/devel   
   module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0
   module load blink-pipeline/main blink_test_data/devel
else
   if [[ -n $module_path ]]; then
      # module purge
      # module load gcc/8.3.0 cascadelake blink_test_data/devel blink-correlator/devel-cmplx blink-imager/cristian-dev
      module reset 
      module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0
      module load blink-pipeline/main blink_test_data/devel
   else
      echo "INFO : non-HPC environment -> modules ignored"
   fi
fi
# salloc --partition gpuq --time 1:00:00 --nodes=1

datapath=$BLINK_TEST_DATADIR/mwa/1276619416/voltages/
if [[ -n "$1" && "$1" != "-" ]]; then
   datapath="$1"
fi

echo "ln -sf ${datapath}/1276619416_1276619418_ch133.dat"
ln -sf ${datapath}/1276619416_1276619418_ch133.dat

if [[ ! -s 20200619163000.metafits ]]; then
   echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits ."
   cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits .
else
   echo "INFO : file 20200619163000.metafits already exists"
fi   

if [[ ! -s calsolutions_chan0_xx.txt ]]; then
   echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan0_xx.txt ."
   cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan0_xx.txt .
else
   echo "INFO : file calsolutions_chan0_xx.txt already exists"
fi

# No averaging for the start:
# Setonix full path to CPU version (if need to be sure): /software/projects/director2183/msok/blink_pipeline/cpu/pipeline/build/blink_pipeline
echo "------------------------------------------"
echo "blink_pipeline -c 4 -C 0 -t 1.00s -o allcorr/ -n 8192 -f 169.60 -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 100 -r -L -G -s calsolutions_chan0_xx.txt -r -V 100 1276619416_1276619418_ch133.dat"
blink_pipeline -c 4 -C 0 -t 1.00s -o allcorr/ -n 8192 -f 169.60 -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 100 -r -L -G -s calsolutions_chan0_xx.txt -r -V 100 1276619416_1276619418_ch133.dat 
echo "------------------------------------------"

