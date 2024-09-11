#!/bin/bash

# TODO : could be sorted by modules/environment if need be:
pipeline_path=/software/projects/director2183/msok/blink_pipeline/gpu/pipeline/

n_avg=32 # 32 fine channels of 40kHz -> 1 coarse channel 
if [[ -n "$1" && "$1" != "-" ]]; then
   n_avg=$1
fi

n_times=10
if [[ -n "$2" && "$2" != "-" ]]; then
   n_times=$2
fi

t=0
while [[ $t -lt $n_times ]];
do
   t_str=`echo $t | awk '{printf("%06d\n",$1);}'`
   ls ch???/test_image_time${t_str}_ch?????_real.fits > fits_list_${t_str}
   
   echo "sbatch ${pipeline_path}/experiments/avg_n_images.sh ${n_avg} fits_list_${t_str} ${t_str}"
   sbatch ${pipeline_path}/experiments/avg_n_images.sh ${n_avg} fits_list_${t_str} ${t_str}
   
   t=$(($t+1))
done

# ls avg_time??????.fits > fits_list_time
# echo "sbatch ${pipeline_path}/experiments/avg_images.sh fits_list_time avg.fits rms.fits"
# sbatch ${pipeline_path}/experiments/avg_images.sh fits_list_time avg.fits rms.fits
