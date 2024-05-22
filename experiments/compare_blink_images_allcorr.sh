#!/bin/bash

refdir=/data_archive/scratch/msok/mwa/1276619416/1276619416_gara/FluxDensity_comparison/casadump/cablecorr_geocorr/000000/
if [[ -n "$1" && "$1" != "-" ]]; then
   refdir="$1"
fi

for chdir in `ls -d ch???`
do
   ch_str=`echo $chdir | cut -b 3-6`
   cd ${chdir}
   
   # cable/geo corrected correlation matrix:
   echo "calcfits_bg xcorr_re.fits = ${refdir}/1276619416_20200619163000_vis_real_channel${ch_str}_time000000_pol0.fits -p 0.5 -p 0.5"
   calcfits_bg xcorr_re.fits = ${refdir}/1276619416_20200619163000_vis_real_channel${ch_str}_time000000_pol0.fits -p 0.5 -p 0.5
   
   echo "calcfits_bg xcorr_im.fits = ${refdir}/1276619416_20200619163000_vis_imag_channel${ch_str}_time000000_pol0.fits -p 0.5 -p 0.5"
   calcfits_bg xcorr_im.fits = ${refdir}/1276619416_20200619163000_vis_imag_channel${ch_str}_time000000_pol0.fits -p 0.5 -p 0.5

   # UV grid counter and real/imag:      
   echo "calcfits_bg uv_grid_counter_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_counter_8192x8192.fits"
   calcfits_bg uv_grid_counter_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_counter_8192x8192.fits

   echo "calcfits_bg uv_grid_real_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_real_8192x8192.fits -p 0.5 -p 0.5"
   calcfits_bg uv_grid_real_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_real_8192x8192.fits -p 0.5 -p 0.5
   
   echo "calcfits_bg uv_grid_imag_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_imag_8192x8192.fits -p 0.5 -p 0.5"
   calcfits_bg uv_grid_imag_8192x8192.fits = ${refdir}/${ch_str}perfreq8192x8192/uv_grid_imag_8192x8192.fits -p 0.5 -p 0.5
   
   # real/imag sky image (FFT 2D):
   echo "calcfits_bg test_image_time000000_ch00${ch_str}_real.fits = ${refdir}/${ch_str}perfreq8192x8192/dirty_image*real.fits -p 0.35 -p 0.35"
   calcfits_bg test_image_time000000_ch00${ch_str}_real.fits = ${refdir}/${ch_str}perfreq8192x8192/dirty_image*real.fits -p 0.35 -p 0.35   
   
   echo "calcfits_bg test_image_time000000_ch00${ch_str}_imag.fits = ${refdir}/${ch_str}perfreq8192x8192/dirty_image*imag.fits -p 0.35 -p 0.35"
   calcfits_bg test_image_time000000_ch00${ch_str}_imag.fits = ${refdir}/${ch_str}perfreq8192x8192/dirty_image*imag.fits -p 0.35 -p 0.35   
      
   cd ..
done
