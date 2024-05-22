#!/bin/bash

refdir=/data_archive/scratch/msok/mwa/1276619416/1276619416_gara/FluxDensity_comparison/casadump/cablecorr_geocorr/000000/
if [[ -n "$1" && "$1" != "-" ]]; then
   refdir="$1"
fi

limit_vis=0.6
if [[ -n "$2" && "$2" != "-" ]]; then
   limit_vis=$2
fi

skyimage_limit=0.005
if [[ -n "$3" && "$3" != "-" ]]; then
   skyimage_limit=$3
fi


for chdir in `ls -d ch???`
do
   ch_str=`echo $chdir | cut -b 3-6`
   cd ${chdir}

   # Because apply solutions puts ZERO on diagonal in correlation matrix (both REAL and IMAG).
   # I use separate script to compare and exclude diagonal:   
   # cable/geo corrected correlation matrix:
   echo "python /home/msok/github/blink/pipeline/experiments/comparefits.py xcorr_re.fits ${refdir}/1276619416_20200619163000_vis_real_channel${ch_str}_time000000_pol0.fits 0 ${limit_vis}"
   python /home/msok/github/blink/pipeline/experiments/comparefits.py xcorr_re.fits ${refdir}/1276619416_20200619163000_vis_real_channel${ch_str}_time000000_pol0.fits 0 ${limit_vis} 

   echo "python /home/msok/github/blink/pipeline/experiments/comparefits.py xcorr_im.fits ${refdir}/1276619416_20200619163000_vis_imag_channel${ch_str}_time000000_pol0.fits 0 ${limit_vis}"
   python /home/msok/github/blink/pipeline/experiments/comparefits.py xcorr_im.fits ${refdir}/1276619416_20200619163000_vis_imag_channel${ch_str}_time000000_pol0.fits 0 ${limit_vis} 
   
   # UV grid counter and real/imag:      
   echo "calcfits_bg uv_grid_counter_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_counter_8192x8192.fits"
   calcfits_bg uv_grid_counter_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_counter_8192x8192.fits

   echo "calcfits_bg uv_grid_real_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_real_8192x8192.fits -p ${limit_vis} -p ${limit_vis}"
   calcfits_bg uv_grid_real_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_real_8192x8192.fits -p ${limit_vis} -p ${limit_vis}
   
   echo "calcfits_bg uv_grid_imag_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_imag_8192x8192.fits -p ${limit_vis} -p ${limit_vis}"
   calcfits_bg uv_grid_imag_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_imag_8192x8192.fits -p ${limit_vis} -p ${limit_vis}
   
   # real/imag sky image (FFT 2D):
   echo "calcfits_bg test_image_time000000_ch00${ch_str}_real.fits = ${refdir}/${ch_str}perfreq8192x8192/dirty_image*real.fits -p ${skyimage_limit} -p ${skyimage_limit}"
   calcfits_bg test_image_time000000_ch00${ch_str}_real.fits = ${refdir}/${ch_str}perfreq8192x8192/dirty_image*real.fits -p ${skyimage_limit} -p ${skyimage_limit}   
   
   echo "calcfits_bg test_image_time000000_ch00${ch_str}_imag.fits = ${refdir}/${ch_str}perfreq8192x8192/dirty_image*imag.fits -p ${skyimage_limit} -p ${skyimage_limit}"
   calcfits_bg test_image_time000000_ch00${ch_str}_imag.fits = ${refdir}/${ch_str}perfreq8192x8192/dirty_image*imag.fits -p ${skyimage_limit} -p ${skyimage_limit}   
      
   cd ..
done
