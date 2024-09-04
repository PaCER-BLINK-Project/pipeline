#!/bin/bash

#SBATCH --account=director2183-gpu # your account pawsey0809-gpu or director2183-gpu
#SBATCH --partition=gpu            # Using the gpu partition
#SBATCH --time=23:59:59
#SBATCH --ntasks-per-node=1        # Set this for 1 mpi task per compute device
#SBATCH --gres=gpu:1
#SBATCH --gpu-bind=closest         # Bind each MPI taks to the nearest GPU
#SBATCH --output=./avg_n_images.o%j
#SBATCH --error=./avg_n_images.e%j
#SBATCH --export=NONE

pipeline_module=blink-pipeline/maingpu
module reset 
module unload gcc/12.2.0 
module swap pawseyenv/2024.05 pawseyenv/2023.08 
module load gcc/12.2.0
module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0 /software/setonix/2023.08/modules/zen3/gcc/12.2.0/libraries
# module load msfitslib/devel 
module load ${pipeline_module} 

n_avg=128
if [[ -n "$1" && "$1" != "-" ]]; then
   n_avg=$1
fi

list=fits_list
if [[ -n "$2" && "$2" != "-" ]]; then
   list=$2
fi

fits_count=`cat $list|wc -l`
start_index=0
while [[ $start_index -lt $fits_count ]];
do
   i_str=`echo $start_index | awk '{printf("%05d\n",$1);}'`
   awk -v n_avg=$n_avg -v start=${start_index} '{if(NR>start && NR<(start+n_avg+1)){print $0;}}' $list > avg_list_${i_str}
   
   echo "avg_images avg_list_${i_str} avg_${i_str}.fits rms_${i_str}.fits -r 10000000.00"
   avg_images avg_list_${i_str} avg_${i_str}.fits rms_${i_str}.fits -r 10000000.00
   
   start_index=$(($start_index+$n_avg))   
done
