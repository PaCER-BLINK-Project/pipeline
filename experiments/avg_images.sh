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

list=fits_list
if [[ -n "$1" && "$1" != "-" ]]; then
   list=$1
fi

avg_image=avg.fits
if [[ -n "$2" && "$2" != "-" ]]; then
   avg_image="$2"
fi

rms_image=rms.fits
if [[ -n "$3" && "$3" != "-" ]]; then
   rms_image="$3"
fi

echo "avg_images $list $avg_image $rms_image -r 10000000.00"
avg_images $list $avg_image $rms_image -r 10000000.00
