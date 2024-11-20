#!/bin/bash

# /home/msok/Desktop/PAWSEY/PaCER/logbook/20230908_testing_BLINK_PIPELINE_cable_and_geom_correction_applycal.odt
# see /media/msok/80f59a5b-8bba-4392-8b31-a4f02fbee2f1/mwa/vcs/1276619416/1276619416_new/cotter/cotter_only_geomandcablecorr_applycal/BLINK_PIPELINE

#SBATCH --account=director2183-gpu # your account pawsey0809-gpu or director2183-gpu
#SBATCH --partition=gpu            # Using the gpu partition
#SBATCH --time=23:59:59
#SBATCH --ntasks-per-node=1        # Set this for 1 mpi task per compute device
#SBATCH --gres=gpu:1
#SBATCH --gpu-bind=closest         # Bind each MPI taks to the nearest GPU
#SBATCH --output=./image_mwa_allch_gpu_setonix.o%j
#SBATCH --error=./image_mwa_allch_gpu_setonix.e%j
#SBATCH --export=NONE

is_one_of () {
    eval "case $1 in ($2) return 0; esac"
    return 1
}


n_channels=768
start_coarse_channel=133
fine_bw=0.04

second=1276619419
if [[ -n "$1" && $1 != "-" ]]; then
   second=$1
fi
ux=$(($second+315964782+1))

obsid=1276619416
if [[ -n "$2" && $2 != "-" ]]; then
   obsid=$2
fi

channel_selection=""
# if [[ -n "$2" && $2 != "-" ]]; then
#   channel_selection="$2"
#fi

if [[ -n "$3" && $3 != "-" ]]; then
   start_coarse_channel=$3
fi

if [[ -s ${obsid}.metafits ]]; then
   echo "INFO : metafits file ${obsid}.metafits found"
else
   echo "WARNING : metafits file ${obsid}.metafits not found -> downloading"
   url="http://ws.mwatelescope.org/metadata/fits?obs_id="
   wget ${url}${obsid} -O ${obsid}.metafits      
fi

# 21,25,58,71,80,81,92,101,108,114,119,125 for obsID = 1276619416
flagged_antennas="19,52,73,78,88,105,107,118,122,125" # for obsID = 1141224136
flag_file=${obsid}_flagged_tiles.txt
if [[ -s ${flag_file} ]]; then
   echo "INFO : found flagged antennas file ${flag_file} -> using these antennas"
   flagged_antennas=`awk '{printf("%s,",$1);}' ${obsid}_flagged_tiles.txt`
else
   echo "WARNING : flag file ${obsid}_flagged_tiles.txt not found - will use default list of antennas to flag $flagged_antennas (for obsID=1141224136)"
fi

if [[ -n "$4" && $4 != "-" ]]; then
   flagged_antennas="$4"
fi

echo "####################################################"
echo "PARAMETERS :"
echo "####################################################"
echo "second = $second ( ux = $ux )"
echo "obsid  = $obsid" 
echo "flagged_antennas = $flagged_antennas"
echo "####################################################"


# pipeline_path=/media/msok/5508b34c-040a-4dce-a8ff-2c4510a5d1a3/pacer/software/blink_pipeline/pipeline/build_20240501_cpu_external/blink_pipeline
# if [[ -n "$1" && "$1" == "setonix" ]] ;then
pipeline_path=/software/projects/director2183/msok/blink_pipeline/gpu/pipeline/build_gpu/blink_pipeline
# fi

module_path=`which module`

if [[ $PAWSEY_CLUSTER = "setonix" ]]; then
   pipeline_module=blink-pipeline/maingpu
   
   module reset 
   module unload gcc/12.2.0 
   module swap pawseyenv/2024.05 pawseyenv/2023.08 
   module load gcc/12.2.0
   
   module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0 /software/setonix/2023.08/modules/zen3/gcc/12.2.0/libraries
   module load msfitslib/devel   
   module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0
   # module load ${pipeline_module} blink_test_data/devel
   
   module use /software/setonix/unsupported
   # module load rocm/5.7.3   
   module load blink_test_data/devel  blink_astroio_msok/master  blink_correlator_msok/gpu blink-imager-gpu/gpu blink_preprocessing/cpu rocm/5.7.3
fi
# salloc --partition gpuq --time 1:00:00 --nodes=1

# echo "cp $BLINK_TEST_DATADIR/mwa/${obsid}/calsolutions_chan*_xx.txt ."
# cp $BLINK_TEST_DATADIR/mwa/${obsid}/calsolutions_chan*_xx.txt .

echo "------------------- preparation completed -------------------"

coarse_channel=$start_coarse_channel
end_coarse_channel=$(($start_coarse_channel+24))

while [[ $coarse_channel -lt $end_coarse_channel ]];
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
 
   echo "DEBUG : checking file /scratch/mwavcs/msok/${obsid}/combined/${obsid}_${second}_ch${coarse_channel}.dat ..."
   if [[ -s /scratch/mwavcs/msok/${obsid}/combined/${obsid}_${second}_ch${coarse_channel}.dat ]]; then   
      echo "ln -sf /scratch/mwavcs/msok/${obsid}/combined/${obsid}_${second}_ch${coarse_channel}.dat"
      ln -sf /scratch/mwavcs/msok/${obsid}/combined/${obsid}_${second}_ch${coarse_channel}.dat
   else
      echo "WARNING : file /scratch/mwavcs/msok/${obsid}/combined/${obsid}_${second}_ch${coarse_channel}.dat not found"
      if [[ -s $BLINK_TEST_DATADIR/mwa/${obsid}/voltages/${obsid}_${second}_ch${coarse_channel}.dat ]]; then
         echo "ln -sf $BLINK_TEST_DATADIR/mwa/${obsid}/voltages/${obsid}_${second}_ch${coarse_channel}.dat"
         ln -sf $BLINK_TEST_DATADIR/mwa/${obsid}/voltages/${obsid}_${second}_ch${coarse_channel}.dat      
      else
         echo "ERROR : missing file $BLINK_TEST_DATADIR/mwa/${obsid}/voltages/${obsid}_${second}_ch${coarse_channel}.dat -> cannot continue"
         exit;
      fi
   fi
   
   if [[ ! -s ${obsid}_${second}_ch${coarse_channel}.dat ]]; then
      echo "ERROR : file ${obsid}_${second}_ch${coarse_channel}.dat does not exist -> cannot continue"
      exit;
   fi

   # this is to pick up correct calibration solution , in the code negative fine channel is interpreted as start fine channel :
   coarse_channel_idx=$(($coarse_channel-$start_coarse_channel))
   start_fine_channel=$(($coarse_channel_idx*32))
   start_fine_channel_param=-1
   if [[ coarse_channel_idx -gt 0 ]]; then
      start_fine_channel_param=-${start_fine_channel}
   fi

   # gps2ux! 1276619419
   # 1592584201
   # was -U 1592584240
   verbose=100
   echo "time $pipeline_path -c 4 -C ${start_fine_channel_param} -t 1.00s -o ch -n 8192 -f ${freq_mhz} -F 30 -M ${obsid}.metafits -u -U ${ux}  -w N -v ${verbose} -r -L -G -s calsolutions -r -V ${verbose} -A ${flagged_antennas} ${obsid}_${second}_ch${coarse_channel}.dat > ${coarse_channel}.out 2>&1"
   time $pipeline_path -c 4 -C ${start_fine_channel_param} -t 1.00s -o ch -n 8192 -f ${freq_mhz} -F 30 -M ${obsid}.metafits -u -U ${ux} -w N -v ${verbose} -r -L -G -s calsolutions -r -V ${verbose} -A ${flagged_antennas} ${obsid}_${second}_ch${coarse_channel}.dat > ${coarse_channel}.out 2>&1
   
   coarse_channel=$(($coarse_channel+1))
done

echo "------------------- blink_pipeline execution completed -------------------"

echo "Create average image"
module load msfitslib/devel
ls ch???/test_image_time000000_ch?????_real.fits > fits_list_all
echo "avg_images fits_list_all avg_all.fits rms_all.fits -r 10000000.00"
avg_images fits_list_all avg_all.fits rms_all.fits -r 10000000.00
