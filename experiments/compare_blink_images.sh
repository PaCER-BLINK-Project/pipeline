#!/bin/bash

#SBATCH --account=director2183
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64gb
#SBATCH --output=./compare_blink_images.o%j
#SBATCH --error=./compare_blink_images.e%j
#SBATCH --export=NONE

module load msfitslib/devel

reference_path=/scratch/mwavcs/msok/1276619416_repeat/Setonix/blink/gpu_casa_dump/nimbus4/
if [[ -n "$1" && "$1" != "-" ]]; then
   reference_path="$1"
fi

postfix="perfreq8192x8192"
if [[ $# -ge 2 && "$2" != "-" ]]; then
   postfix="$2"
fi

prefix=""
if [[ $# -ge 3 && "$3" != "-" ]]; then
   prefix="$3"
fi

for chdir in `ls -d ${prefix}???${postfix}`
do
   echo "Comparing data in ${chdir} vs. ${reference_path}/${chdir}/ :"
   if [[ -d ${reference_path}/${chdir} ]]; then
      cd ${chdir}/
      for fits in `ls *.fits`; do    
#          echo "calcfits_bg $fits = ${reference_path}/${fits} | grep \"COMPARISON RESULT\""
          # dirty_image_20240404T081559000_real.fits
          fits_ref="${reference_path}/${chdir}/${fits}"
          fits_ref_prefix=`echo $fits | cut -b 1-11`
          if [[ $fits_ref_prefix == "dirty_image" ]]; then
             fits_ref_postfix=`echo $fits | cut -b 31-`
             echo "ls ${reference_path}/${chdir}/${fits_ref_prefix}*${fits_ref_postfix}"
             fits_ref=`ls ${reference_path}/${chdir}/${fits_ref_prefix}*${fits_ref_postfix} | head -1`          
          fi
          
          if [[ $fits == "vis_re_cable_corr_gpu.fits" ]]; then
             if [[ ! -s ${fits_ref} ]]; then
                fits_ref=${reference_path}/${chdir}/vis_re_cable_corr.fits
             fi
          fi

          if [[ $fits == "vis_im_cable_corr_gpu.fits" ]]; then
             if [[ ! -s ${fits_ref} ]]; then
                fits_ref=${reference_path}/${chdir}/vis_im_cable_corr.fits
             fi
          fi

          echo "calcfits_bg $fits = ${fits_ref}"
#          calcfits_bg $fits = ${fits_ref} | grep "COMPARISON RESULT for file $fits"
          calcfits_bg $fits = ${fits_ref} | grep "COMPARISON RESULT" 
          echo
          sleep 5
      done
      echo "----------------------------------------"
      echo
      echo
      sleep 2;      
      cd -
   else
      echo "Skipped : ${reference_path}/${chdir}"
   fi
done
