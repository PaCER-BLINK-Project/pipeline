#!/bin/bash

n_seconds=1
if [[ -n "$1" && "$1" != "-" ]]; then
   n_seconds=$1
fi

obsid=1276619416
if [[ -n "$2" && "$2" != "-" ]]; then
   obsid=$2
fi

start_second=1276619418
if [[ -n "$3" && "$3" != "-" ]]; then
   start_second=$3
fi

avg_channels=4
if [[ -n "$4" && $4 != "-" ]]; then
   avg_channels=$4
fi


url="http://ws.mwatelescope.org/metadata/fits?obs_id="
pipeline_path=/software/projects/director2183/msok/blink_pipeline/gpu/pipeline/

second=$start_second
end_second=$(($start_second+$n_seconds))

while [[ $second -lt $end_second ]];
do
   subdir=time_${second}
   mkdir -p ${subdir}
    
   cd ${subdir}
   
   if [[ -s ${obsid}.metafits ]]; then
      echo "Metafits ${obsid}.metafits found OK"
   else
      echo "wget ${url}${obsid} -O ${obsid}.metafits"
      wget ${url}${obsid} -O ${obsid}.metafits
      
      if [[ ! -s ${obsid}.metafits ]]; then
         echo "cp ../template/${obsid}.metafits ."
         cp ../template/${obsid}.metafits .
      fi
   fi
   
   echo "sbatch ${pipeline_path}/experiments/test_mwa_obsid1276619416_full_imaging_loop_chall_autometafits.sh - 100 ${avg_channels} $obsid $second"
   sbatch ${pipeline_path}/experiments/test_mwa_obsid1276619416_full_imaging_loop_chall_autometafits.sh - 100 ${avg_channels} $obsid $second
   cd ..
   
   second=$(($second+1))
done
   
   
   
