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

reference_path=/data_archive/scratch/msok/mwa/1276619416/1276619416_gara/FluxDensity_comparison/blink/allcorr/setonix/
if [[ -n "$1" && "$1" != "-" ]]; then
   reference_path="$1"
fi

postfix=""
if [[ $# -ge 2 && "$2" != "-" ]]; then
   postfix="$2"
fi

prefix="ch"
if [[ $# -ge 3 && "$3" != "-" ]]; then
   prefix="$3"
fi

refpostfix=""
gzipped=0
if [[ $# -ge 4 && "$4" != "-" ]]; then
   gzipped=$4   
fi
if [[ $gzipped -gt 0 ]]; then
   refpostfix=".gz"
fi

compare_avg=1
if [[ $# -ge 5 && "$5" != "-" ]]; then
   compare_avg=$5
fi


for chdir in `ls -d ${prefix}???${postfix}`
do
   echo "Comparing data in ${chdir} vs. ${reference_path}/${chdir}/ :"
   if [[ -d ${reference_path}/${chdir} ]]; then
      cd ${chdir}/
      refdir=${reference_path}/${chdir}/
      
      if [[ -d ${refdir} ]]; then
         for fits in `ls *.fits`; do    
#             echo "calcfits_bg $fits = ${reference_path}/${fits} | grep \"COMPARISON RESULT\""
             # dirty_image_20240404T081559000_real.fits
             fits_ref="${reference_path}/${chdir}/${fits}${refpostfix}"
             fits_ref_prefix=`echo $fits | cut -b 1-11`
          
             if [[ -s ${fits_ref} ]]; then
                echo "calcfits_bg $fits = ${fits_ref} -p 0.00005 -p 0.00005 -X"
                calcfits_bg $fits = ${fits_ref} -p 0.00005 -p 0.00005 -X 
             else
                echo "WARNING : file ${fits_ref} does not exist (reference file for ${fits} not found) -> comparison skipped"
             fi
             echo
             sleep 5
         done
         echo "----------------------------------------"
         echo
         echo
      else
         echo "INFO : reference directory ${refdir} does not exist -> skipped"
      fi
      cd -
   else
      echo "Skipped : ${reference_path}/${chdir}"
   fi
done

if [[ $compare_avg -gt 0 ]]; then
   echo "calcfits_bg avg_all.fits = ${reference_path}/avg_all.fits"
   calcfits_bg avg_all.fits = ${reference_path}/avg_all.fits
fi   

