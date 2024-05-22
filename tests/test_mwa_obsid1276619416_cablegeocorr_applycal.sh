#!/bin/bash

# /home/msok/Desktop/PAWSEY/PaCER/logbook/20230908_testing_BLINK_PIPELINE_cable_and_geom_correction_applycal.odt
# see /media/msok/80f59a5b-8bba-4392-8b31-a4f02fbee2f1/mwa/vcs/1276619416/1276619416_new/cotter/cotter_only_geomandcablecorr_applycal/BLINK_PIPELINE

module_path=`which module`

if [[ $PAWSEY_CLUSTER = "setonix" ]]; then
   module reset
   module load blink_test_data/devel  blink_astroio/master  blink_correlator/master blink-imager-cpu/devel blink_preprocessing/master
else
   if [[ -n $module_path ]]; then
      module purge
      module load gcc/8.3.0 cascadelake blink_test_data/devel blink-correlator/devel-cmplx blink-imager/cristian-dev
   else
      echo "INFO : non-HPC environment -> modules ignored"
   fi
fi
# salloc --partition gpuq --time 1:00:00 --nodes=1

echo "ln -sf $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch133.dat"
ln -sf $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch133.dat

if [[ ! -s 20200619163000.metafits ]]; then
   echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits ."
   cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits .
else
   echo "INFO : file 20200619163000.metafits already exists"
fi   

if [[ ! -s calsolutions_chan0_xx.txt ]]; then
   echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan0_xx.txt ."
   cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan0_xx.txt .
else
   echo "INFO : file calsolutions_chan0_xx.txt already exists"
fi

echo "------------------------------------------"
echo "./blink_pipeline -c 4 -C 0 -t 1.00s -o allcorr/ -n 2048 -f 169.60 -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 100 -r -L -G -s calsolutions_chan0_xx.txt -r -V 100 1276619416_1276619418_ch133.dat"
./blink_pipeline -c 4 -C 0 -t 1.00s -o allcorr/ -n 2048 -f 169.60 -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 100 -r -L -G -s calsolutions_chan0_xx.txt -r -V 100 1276619416_1276619418_ch133.dat 
echo "------------------------------------------"

max_diff=0.001 # 0.4/0.3 to compare with cotter files 1276619416_20200619163000_vis_real_channel000_time000000_pol0.fits
echo "calcfits_bg cal_re.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/cal_re.fits -p ${max_diff} -p ${max_diff}"
calcfits_bg cal_re.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/cal_re.fits -p ${max_diff} -p ${max_diff}
# repeat to get exit code :
exit_code=0
equal=`calcfits_bg cal_re.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/cal_re.fits -p ${max_diff} -p ${max_diff} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi


echo "calcfits_bg cal_im.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/cal_im.fits -p ${max_diff} -p ${max_diff}"
calcfits_bg cal_im.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/cal_im.fits -p ${max_diff} -p ${max_diff}
# repeat to get exit code :
exit_code=0
equal=`calcfits_bg cal_im.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/cal_im.fits -p ${max_diff} -p ${max_diff} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

# compare images:
max_diff=0.00001
echo "calcfits_bg allcorr/test_image_time000000_ch00000_real.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/test_image_real.fits -p ${max_diff} -p ${max_diff}"
calcfits_bg allcorr/test_image_time000000_ch00000_real.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/test_image_real.fits -p ${max_diff} -p ${max_diff}
exit_code=0
equal=`calcfits_bg allcorr/test_image_time000000_ch00000_real.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/test_image_real.fits -p ${max_diff} -p ${max_diff} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

echo "calcfits_bg allcorr/test_image_time000000_ch00000_imag.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/test_image_imag.fits -p ${max_diff} -p ${max_diff}"
calcfits_bg allcorr/test_image_time000000_ch00000_imag.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/test_image_imag.fits -p ${max_diff} -p ${max_diff}
exit_code=0
equal=`calcfits_bg allcorr/test_image_time000000_ch00000_imag.fits = $BLINK_TEST_DATADIR/mwa/1276619416/results/cablegeocorr_applycal/test_image_imag.fits -p ${max_diff} -p ${max_diff} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

me=$(basename "$0")
echo "Exiting script $me with exit code = $exit_code"
exit $exit_code

# Comparisons with cotter calibrated:
# calcfits_bg cal_re.fits = ../1276619416_20200619163000_vis_real_channel000_time000000_pol0.fits -p 0.1 -p 0.1
# COMPARISON RESULT : Images differ by 704 pixels
# calcfits_bg cal_im.fits = ../1276619416_20200619163000_vis_imag_channel000_time000000_pol0.fits -p 0.1 -p 0.1
# COMPARISON RESULT : Images differ by 580 pixels
