#!/bin/bash

# see flowchart - for what to compare with what !!!
# /home/msok/Desktop/PAWSEY/PaCER/logbook/BLINK_pipeline_xcorr_version.odt
# /home/msok/Desktop/EDA2/papers/2024/BLINK_pipeline/images/BLINK_pipeline_xcorr_version.odg

#SBATCH --account=director2183
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64gb
#SBATCH --output=./compare_blink_images_gpu_vs_cpu.o%j
#SBATCH --error=./compare_blink_images_gpu_vs_cpu.e%j
#SBATCH --export=NONE

module load msfitslib/devel

reference_path=/scratch/mwavcs/msok/1276619416_repeat/Setonix/blink/gpu_casa_dump/nimbus4/
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

for chdir in `ls -d ${prefix}???${postfix}`
do
   echo "Comparing data in ${chdir} vs. ${reference_path}/${chdir}/ :"
   if [[ -d ${reference_path}/${chdir} ]]; then
      cd ${chdir}/
      # manual process because not always the same name FITS file corresponds to the same stage of processing in GPU vs. CPU versions
      # For example, xcorr_re.fits in CPU version is after applying all: cable, geo. corr and cal. solutions
      #              but in GPU versions it's only after applying cal. solutions as cable and geo. corrections are applied later in GPU kernels in gridding_fast(...) :
      echo "---------------------------------------- $chdir ----------------------------------------"
      # self consistency of GPU version:
      echo "calcfits_bg re.fits = test_vis_re.fits"
      calcfits_bg re.fits = test_vis_re.fits      
      
      echo "calcfits_bg im.fits = test_vis_im.fits"
      calcfits_bg im.fits = test_vis_im.fits

      # GPU vs. CPU :
      # gridding_imaging_gpu_precorr - GPU version pre cable/geo corr. after calsol. only, no counterpart in CPU version -> nothing to compare to.
      
      # after correlation files:
      echo "calcfits_bg re.fits = ${reference_path}/${chdir}/re.fits"
      calcfits_bg re.fits = ${reference_path}/${chdir}/re.fits
      
      echo "calcfits_bg test_vis_re.fits = ${reference_path}/${chdir}/test_vis_re.fits"
      calcfits_bg test_vis_re.fits = ${reference_path}/${chdir}/test_vis_re.fits
            
      # UVW files (u,v,w).fits:      
      for uvw_fits in `ls ?.fits`
      do
         echo "calcfits_bg $uvw_fits = ${reference_path}/${chdir}/${uvw_fits}"
         calcfits_bg $uvw_fits = ${reference_path}/${chdir}/${uvw_fits}
      done

      # UV grids and counter separately as there are different tolarances for EQUALITY :      
      echo "calcfits_bg uv_grid_real_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_real_8192x8192.fits -p 0.5 -p 0.5"
      calcfits_bg uv_grid_real_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_real_8192x8192.fits -p 0.5 -p 0.5
      
      echo "calcfits_bg uv_grid_imag_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_imag_8192x8192.fits -p 0.5 -p 0.5"
      calcfits_bg uv_grid_imag_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_imag_8192x8192.fits -p 0.5 -p 0.5
      
      echo "calcfits_bg uv_grid_counter_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_counter_8192x8192.fits"
      calcfits_bg uv_grid_counter_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_counter_8192x8192.fits
      
      # GPU after cable correction is the same as xcorr_re/im in CPU version:
      for fits in `ls vis_??_cable_corr_gpu.fits`
      do
         echo "calcfits_bg ${fits} = ${reference_path}/${chdir}/${fits}"
         calcfits_bg ${fits} = ${reference_path}/${chdir}/${fits}
      done
      
      # Geo corr:      
      for fits in `ls vis_??_geom_corr.fits`
      do
         echo "calcfits_bg ${fits} = ${reference_path}/${chdir}/${fits}"
         calcfits_bg ${fits} = ${reference_path}/${chdir}/${fits}
      done

      # xcorr :
      for fits in `ls xcorr_??.fits`
      do
         echo "calcfits_bg ${fits} = ${reference_path}/${chdir}/${fits}"
         calcfits_bg ${fits} = ${reference_path}/${chdir}/${fits}
      done
      
      for fits in `ls dirty_test_????_8192x8192.fits`
      do
         echo "calcfits_bg ${fits} = ${reference_path}/${chdir}/${fits}"
         calcfits_bg ${fits} = ${reference_path}/${chdir}/${fits}
      done
      
      # final 2D FFT (sky images):
      sky_real=`ls test_image_time??????_ch?????_real.fits  | tail -1`
      echo "calcfits_bg ${sky_real} = ${reference_path}/${chdir}/${sky_real} -p 0.004 -p 0.004"
      calcfits_bg ${sky_real} = ${reference_path}/${chdir}/${sky_real} -p 0.004 -p 0.004
      
      sky_imag=`ls test_image_time??????_ch?????_imag.fits  | tail -1`
      echo "calcfits_bg ${sky_imag} = ${reference_path}/${chdir}/${sky_imag} -p 0.00001 -p 0.00001"
      calcfits_bg ${sky_imag} = ${reference_path}/${chdir}/${sky_imag} -p 0.00001 -p 0.00001
      
      echo "----------------------------------------------------------------------------------------"
      echo
      echo
      sleep 2;      
      cd -
   else
      echo "Skipped : ${reference_path}/${chdir}"
   fi
done
