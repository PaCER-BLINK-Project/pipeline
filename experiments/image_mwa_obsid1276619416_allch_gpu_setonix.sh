#!/bin/bash

#SBATCH --account=director2183-gpu # your account pawsey0809-gpu or director2183-gpu
#SBATCH --partition=gpu            # Using the gpu partition
#SBATCH --time=06:00:00
#SBATCH --ntasks-per-node=1        # Set this for 1 mpi task per compute device
#SBATCH --gres=gpu:1
#SBATCH --gpu-bind=closest         # Bind each MPI taks to the nearest GPU
#SBATCH --output=./image_mwa_obsid1276619416_allch_gpu_setonix.o%j
#SBATCH --error=./image_mwa_obsid1276619416_allch_gpu_setonix.e%j
#SBATCH --export=NONE

# /home/msok/Desktop/PAWSEY/PaCER/logbook/20230908_testing_BLINK_PIPELINE_cable_and_geom_correction_applycal.odt
# see /media/msok/80f59a5b-8bba-4392-8b31-a4f02fbee2f1/mwa/vcs/1276619416/1276619416_new/cotter/cotter_only_geomandcablecorr_applycal/BLINK_PIPELINE

n_channels=768
start_coarse_channel=133
fine_bw=0.04
n_channels=768

# blink-pipeline/main_ongpu - ongpu 
pipeline_module=blink-pipeline/maingpu
if [[ -n "$4" && "$4" != "-" ]]; then
   pipeline_module="$4"
fi

if [[ $PAWSEY_CLUSTER = "setonix" ]]; then
   module reset 
   module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0 /software/setonix/2023.08/modules/zen3/gcc/12.2.0/libraries
   module load msfitslib/devel   
   module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0
   module load ${pipeline_module} blink_test_data/devel
else
   if [[ -n $module_path ]]; then
      module purge
      module load gcc/8.3.0 cascadelake blink_test_data/devel blink-correlator/devel-cmplx blink-imager/cristian-dev
   else
      echo "INFO : non-HPC environment -> modules ignored"
   fi
fi
# salloc --partition gpuq --time 1:00:00 --nodes=1

datapath=/scratch/mwavcs/msok/1276619416/combined/
if [[ -n "$1" && "$1" != "-" ]]; then
   datapath="$1"
fi

start_coarse_channel=133
if [[ -n "$2" && "$2" != "-" ]]; then
   start_coarse_channel=$2
fi

cotter_compatible=0
if [[ -n "$3" && "$3" != "-" ]]; then
   cotter_compatible=$3
fi

if [[ -n "$5" && "$5" != "-" ]]; then
   n_channels=$5
fi

echo "###########################################"
echo "PARAMETERS:"
echo "###########################################"
echo "datapath = $datapath"
echo "start_coarse_channel = $start_coarse_channel"
echo "cotter_compatible = $cotter_compatible"
echo "pipeline_module   = $pipeline_module"
echo "n_channels = $n_channels"
echo "###########################################"


if [[ ! -s 20200619163000.metafits ]]; then
   echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits ."
   cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits .
else
   echo "INFO : file 20200619163000.metafits already exists"
fi   

# if [[ ! -s calsolutions_chan0_xx.txt ]]; then
#   echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan0_xx.txt ."
#   cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan0_xx.txt .
#else
#   echo "INFO : file calsolutions_chan0_xx.txt already exists"
#fi

# cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan???_xx.txt .
# cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan???_yy.txt .

# No averaging for the start:
echo "------------------------------------------"
ch=0
while [[ $ch -lt $n_channels ]]; 
do
   ch_str=`echo $ch | awk '{printf("%03d\n",$1);}'`
   coarse_channel=`echo $ch | awk -v cc=${start_coarse_channel} '{printf("%d\n",cc+($1/32));}'`
   fine_channel_to_image=`echo $ch | awk -v cc=${start_coarse_channel} '{printf("%d\n",($1 % 32));}'`
   
   echo "Processing fine_channel = $ch -> coarse_channel = $coarse_channel -> fine channel to image = $fine_channel_to_image"
   # WARNING to make it the same as cotter :
   # cotter must have options : -norfi -nosbgains -sbpassband /media/msok/5508b34c-040a-4dce-a8ff-2c4510a5d1a3/mwa/data/1276619416/cotter_laptop/no_cable/subband-passband-32ch-unitary.txt
   # and freq_mhz has to be calculated without part 0.04/2.00
   if [[ $cotter_compatible -gt 0 ]]; then
      freq_mhz=`echo $ch | awk -v start_coarse_channel=${start_coarse_channel} -v fine_bw=${fine_bw} '{freq_mhz=start_coarse_channel*1.28-0.64+fine_bw*$1;printf("%.4f\n",freq_mhz);}'`
   else
      freq_mhz=`echo $ch | awk -v start_coarse_channel=${start_coarse_channel} -v fine_bw=${fine_bw} '{freq_mhz=start_coarse_channel*1.28-0.64+fine_bw/2.00+fine_bw*$1;printf("%.4f\n",freq_mhz);}'`
   fi
   
   echo "ln -sf ${datapath}/1276619416_1276619418_ch${coarse_channel}.dat"
   ln -sf ${datapath}/1276619416_1276619418_ch${coarse_channel}.dat


   mkdir -p ch${ch_str}/ 
   echo "Imaging channel = $ch ($ch_str)"
   # !!!! WARNING : do not touch the explicit path it is here to make sure I am really using the version I want !!!!
   #                trespassers will be punished !!!   
   echo "/software/projects/director2183/msok/blink_pipeline/gpu/pipeline/build_gpu/blink_pipeline -c 4 -C $fine_channel_to_image -t 1.00s -o ch${ch_str}/ -n 8192 -f ${freq_mhz} -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 100 -r -L -G -s calsolutions_chan${ch_str}_xx.txt -r -V 100 1276619416_1276619418_ch${coarse_channel}.dat > ch${ch_str}/out 2>&1"
   /software/projects/director2183/msok/blink_pipeline/gpu/pipeline/build_gpu/blink_pipeline -c 4 -C $fine_channel_to_image -t 1.00s -o ch${ch_str}/ -n 8192 -f ${freq_mhz} -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 100 -r -L -G -s calsolutions_chan${ch_str}_xx.txt -r -V 100 1276619416_1276619418_ch${coarse_channel}.dat > ch${ch_str}/out 2>&1

   # moving test correlation matrix file to specific channel subdirectory:   
   mv ??.fits ch${ch_str}/
   mv *test_vis_??.fits ch${ch_str}/
   mv xcorr_??.fits ch${ch_str}/
   mv vis_??_cable_corr.fits ch${ch_str}/
   mv vis_??_geom_corr.fits ch${ch_str}/
   mv vis_??_cable_corr.fits ch${ch_str}/
   mv gridding_imaging_gpu_precorr_??.fits ch${ch_str}/
   mv vis_??_cable_corr_gpu.fits ch${ch_str}/
   mv vis_??_geom_corr_gpu.fits ch${ch_str}/
   mv corrmatrix_after_cable_and_geo_corr_gpu_??.fits ch${ch_str}/
   mv corrmatrix_before_run_imager_nocorrections_gpu_re.fits ch${ch_str}/
   

   ch=$(($ch+1))
done

echo "Create average image"
module load msfitslib/devel
ls ch???/test_image_time000000_ch?????_real.fits > fits_list_all
echo "avg_images fits_list_all avg_all.fits rms_all.fits -r 10000000.00"
avg_images fits_list_all avg_all.fits rms_all.fits -r 10000000.00
   
echo "------------------------------------------"

