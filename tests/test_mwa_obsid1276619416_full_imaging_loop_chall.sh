#!/bin/bash

# /home/msok/Desktop/PAWSEY/PaCER/logbook/20230908_testing_BLINK_PIPELINE_cable_and_geom_correction_applycal.odt
# see /media/msok/80f59a5b-8bba-4392-8b31-a4f02fbee2f1/mwa/vcs/1276619416/1276619416_new/cotter/cotter_only_geomandcablecorr_applycal/BLINK_PIPELINE

#SBATCH --account=director2183-gpu # your account pawsey0809-gpu or director2183-gpu
#SBATCH --partition=gpu            # Using the gpu partition
#SBATCH --time=06:00:00
#SBATCH --ntasks-per-node=1        # Set this for 1 mpi task per compute device
#SBATCH --gres=gpu:1
#SBATCH --gpu-bind=closest         # Bind each MPI taks to the nearest GPU
#SBATCH --output=./image_mwa_obsid1276619416_allch_gpu_setonix.o%j
#SBATCH --error=./image_mwa_obsid1276619416_allch_gpu_setonix.e%j
#SBATCH --export=NONE

is_one_of () {
    eval "case $1 in ($2) return 0; esac"
    return 1
}


channel_selection=""
if [[ -n "$1" && $1 != "-" ]]; then
   channel_selection="$1"
fi

n_channels=768
start_coarse_channel=133
fine_bw=0.04

pipeline_path=/software/projects/director2183/msok/blink_pipeline/cpu/pipeline/build/blink_pipeline
if [[ -n "$1" && "$1" == "setonix" ]] ;then
   pipeline_path=/software/projects/director2183/msok/blink_pipeline/gpu/pipeline/build_gpu/blink_pipeline
fi

module_path=`which module`

if [[ $PAWSEY_CLUSTER = "setonix" ]]; then
   # pipeline_module=blink-pipeline/maingpu
   pipeline_module=blink-pipeline/maincpu.lua
   
   module reset 
   
   module unload gcc/12.2.0 
   module swap pawseyenv/2024.05 pawseyenv/2023.08 
   module load gcc/12.2.0
        
   module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0 /software/setonix/2023.08/modules/zen3/gcc/12.2.0/libraries
   module load msfitslib/devel   
   module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0
   module load ${pipeline_module} blink_test_data/devel
fi
# salloc --partition gpuq --time 1:00:00 --nodes=1

if [[ ! -s 20200619163000.metafits ]]; then
   echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits ."
   cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits .
else
   echo "INFO : file 20200619163000.metafits already exists"
fi   

#ch=0
#while [[ $ch -lt 32 ]];
#do
#   ch_str=`echo $ch | awk '{printf("%03d",$1);}'`
#   echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan${ch_str}_xx.txt ."
#   cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan${ch_str}_xx.txt .
#   
#   mkdir -p ch${ch_str}
#   
#   ch=$(($ch+1))
#done   
echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan*_xx.txt ."
cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan*_xx.txt .

echo "------------------- preparation completed -------------------"

coarse_channel=$start_coarse_channel
end_coarse_channel=$(($start_coarse_channel+24))

while [[ $coarse_channel -le $end_coarse_channel ]];
do
   if [[ -n "$channel_selection" ]]; then
      if is_one_of "$coarse_channel" '156'; then
         echo "Processing channel $coarse_channel"
      else
         echo "INFO : skipping channel $coarse_channel - not on list to process ($channel_selection)"
         coarse_channel=$(($coarse_channel+1))
         continue;
      fi
   fi
   freq_mhz=`echo $ch | awk -v coarse_channel=${coarse_channel} -v fine_bw=${fine_bw} '{freq_mhz=coarse_channel*1.28;printf("%.4f\n",freq_mhz);}'`
   echo "INFO : processing coarse_channel = $coarse_channel -> freq_mhz = $freq_mhz"      

   if [[ -s /scratch/mwavcs/msok/1276619416/combined/1276619416_1276619418_ch${coarse_channel}.dat ]]; then   
      echo "ln -sf /scratch/mwavcs/msok/1276619416/combined/1276619416_1276619418_ch${coarse_channel}.dat"
      ln -sf /scratch/mwavcs/msok/1276619416/combined/1276619416_1276619418_ch${coarse_channel}.dat
   else
      if [[ -s $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch${coarse_channel}.dat ]]; then
         echo "ln -sf $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch${coarse_channel}.dat"
         ln -sf $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch${coarse_channel}.dat      
      else
         echo "ERROR : missing file $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch${coarse_channel}.dat -> cannot continue"
         exit;
      fi
   fi
   
   if [[ ! -s 1276619416_1276619418_ch${coarse_channel}.dat ]]; then
      echo "ERROR : file 1276619416_1276619418_ch${coarse_channel}.dat does not exist -> cannot continue"
      exit;
   fi

   # this is to pick up correct calibration solution , in the code negative fine channel is interpreted as start fine channel :
   coarse_channel_idx=$(($coarse_channel-$start_coarse_channel))
   start_fine_channel=$(($coarse_channel_idx*32))
   start_fine_channel_param=-1
   if [[ coarse_channel_idx -gt 0 ]]; then
      start_fine_channel_param=-${start_fine_channel}
   fi
   echo "$pipeline_path -c 4 -C ${start_fine_channel_param} -t 1.00s -o ch -n 8192 -f ${freq_mhz} -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 100 -r -L -G -s calsolutions -r -V 100 -A 21,25,58,71,80,81,92,101,108,114,119,125 1276619416_1276619418_ch${coarse_channel}.dat > ${coarse_channel}.out 2>&1"
   $pipeline_path -c 4 -C ${start_fine_channel_param} -t 1.00s -o ch -n 8192 -f ${freq_mhz} -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 100 -r -L -G -s calsolutions -r -V 100 -A 21,25,58,71,80,81,92,101,108,114,119,125 1276619416_1276619418_ch${coarse_channel}.dat > ${coarse_channel}.out 2>&1
   
   coarse_channel=$(($coarse_channel+1))
done

echo "------------------- blink_pipeline execution completed -------------------"

echo "Create average image"
module load msfitslib/devel
ls ch???/test_image_time000000_ch?????_real.fits > fits_list_all
echo "avg_images fits_list_all avg_all.fits rms_all.fits -r 10000000.00"
avg_images fits_list_all avg_all.fits rms_all.fits -r 10000000.00
