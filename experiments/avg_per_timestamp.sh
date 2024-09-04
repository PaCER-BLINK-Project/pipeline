#!/bin/bash

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
   
   echo "sbatch ./avg_n_images.sh ${n_avg} fits_list_${t_str}"
   sbatch ./avg_n_images.sh ${n_avg} fits_list_${t_str}
   
   t=$(($t+1))
done
