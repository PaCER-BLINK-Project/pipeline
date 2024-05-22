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

postfix="perfreq8192x8192"
if [[ $# -ge 2 && "$2" != "-" ]]; then
   postfix="$2"
fi

prefix=""
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
      # manual process because not always the same name FITS file corresponds to the same stage of processing in GPU vs. CPU versions
      # For example, xcorr_re.fits in CPU version is after applying all: cable, geo. corr and cal. solutions
      #              but in GPU versions it's only after applying cal. solutions as cable and geo. corrections are applied later in GPU kernels in gridding_fast(...) :
      echo "---------------------------------------- $chdir ----------------------------------------"
      # self consistency of GPU version:
      echo "calcfits_bg re.fits = test_vis_re.fits${refpostfix}"
      calcfits_bg re.fits = test_vis_re.fits${refpostfix}
      
      echo "calcfits_bg im.fits = test_vis_im.fits${refpostfix}"
      calcfits_bg im.fits = test_vis_im.fits${refpostfix}

      # GPU vs. CPU :
      # gridding_imaging_gpu_precorr - GPU version pre cable/geo corr. after calsol. only, no counterpart in CPU version -> nothing to compare to.
      
      # after correlation files:
      echo "calcfits_bg re.fits = ${reference_path}/${chdir}/re.fits${refpostfix}"
      calcfits_bg re.fits = ${reference_path}/${chdir}/re.fits${refpostfix}
      
      echo "calcfits_bg test_vis_re.fits = ${reference_path}/${chdir}/test_vis_re.fits${refpostfix}"
      calcfits_bg test_vis_re.fits = ${reference_path}/${chdir}/test_vis_re.fits${refpostfix}
            
      # UVW files (u,v,w).fits:      
      for uvw_fits in `ls ?.fits`
      do
         echo "calcfits_bg $uvw_fits = ${reference_path}/${chdir}/${uvw_fits}${refpostfix}"
         calcfits_bg $uvw_fits = ${reference_path}/${chdir}/${uvw_fits}${refpostfix}
      done

      # UV grids and counter separately as there are different tolarances for EQUALITY :      
      echo "calcfits_bg uv_grid_real_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_real_8192x8192.fits${refpostfix} -p 0.9 -p 0.9"
      calcfits_bg uv_grid_real_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_real_8192x8192.fits${refpostfix} -p 0.9 -p 0.9
      
      echo "calcfits_bg uv_grid_imag_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_imag_8192x8192.fits${refpostfix} -p 0.9 -p 0.9"
      calcfits_bg uv_grid_imag_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_imag_8192x8192.fits${refpostfix} -p 0.9 -p 0.9
      
      echo "calcfits_bg uv_grid_counter_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_counter_8192x8192.fits${refpostfix}"
      calcfits_bg uv_grid_counter_8192x8192.fits = ${reference_path}/${chdir}/uv_grid_counter_8192x8192.fits${refpostfix}
      
      # GPU after cable correction is the same as xcorr_re/im in CPU version:
      echo "calcfits_bg vis_re_cable_corr_gpu.fits = ${reference_path}/${chdir}/xcorr_re.fits${refpostfix} -p 0.7 -p 0.7"
      calcfits_bg vis_re_cable_corr_gpu.fits = ${reference_path}/${chdir}/xcorr_re.fits${refpostfix} -p 0.7 -p 0.7
      
      echo "calcfits_bg vis_im_cable_corr_gpu.fits = ${reference_path}/${chdir}/xcorr_im.fits${refpostfix} -p 0.7 -p 0.7"
      calcfits_bg vis_im_cable_corr_gpu.fits = ${reference_path}/${chdir}/xcorr_im.fits${refpostfix} -p 0.7 -p 0.7
      
      # Geo corr:      
      # There is no counterpart of this file in CPU version.
      # to compare full calibrated data compare vis_re_cable_corr_gpu.fits with xcorr_re.fits (and _im) - which is already done
      # echo "calcfits_bg vis_re_geom_corr_gpu.fits = ${reference_path}/${chdir}/vis_re_geom_corr_gpu.fits"
      # calcfits_bg vis_re_geom_corr_gpu.fits = ${reference_path}/${chdir}/vis_re_geom_corr_gpu.fits
      # echo "calcfits_bg vis_im_geom_corr_gpu.fits = ${reference_path}/${chdir}/vis_im_geom_corr_gpu.fits"
      # calcfits_bg vis_im_geom_corr_gpu.fits = ${reference_path}/${chdir}/vis_im_geom_corr_gpu.fits
      
      # final 2D FFT (sky images):
      sky_real=`ls test_image_time??????_ch?????_real.fits${refpostfix}  | tail -1`
      echo "calcfits_bg ${sky_real} = ${reference_path}/${chdir}/${sky_real} -p 0.006 -p 0.006"
      calcfits_bg ${sky_real} = ${reference_path}/${chdir}/${sky_real} -p 0.006 -p 0.006
      
      sky_imag=`ls test_image_time??????_ch?????_imag.fits${refpostfix}  | tail -1`
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

if [[ $compare_avg -gt 0 ]]; then
   echo "calcfits_bg avg_all.fits = ${reference_path}/avg_all.fits"
   calcfits_bg avg_all.fits = ${reference_path}/avg_all.fits
fi   
