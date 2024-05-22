#!/bin/bash
if [ $PAWSEY_CLUSTER = "setonix" ]; then
module reset
module load blink_test_data/devel blink_preprocessing/master  blink_astroio/master blink_correlator/master blink-imager-cpu/devel
 
else
module purge
module load gcc/8.3.0 cascadelake blink_test_data/devel blink-correlator/devel-cmplx blink-imager/cristian-dev
fi
# salloc --partition gpuq --time 1:00:00 --nodes=1

echo "./blink_pipeline -t 0.28311552s -f159.375 -F180.00 -a ${BLINK_TEST_DATADIR}/eda2/antenna_locations.txt -i eda2 -Z ${BLINK_TEST_DATADIR}/eda2/channel_cont_20220118_41581_0_binary.bin"
./blink_pipeline -t 0.28311552s -f159.375 -F180.00 -a ${BLINK_TEST_DATADIR}/eda2/antenna_locations.txt -i eda2 -Z ${BLINK_TEST_DATADIR}/eda2/channel_cont_20220118_41581_0_binary.bin


echo "calcfits_bg test_image_real.fits = $BLINK_TEST_DATADIR/eda2/test_image_real.fits"
calcfits_bg test_image_real.fits = $BLINK_TEST_DATADIR/eda2/test_image_real.fits
# repeat to get exit code :
exit_code=0
equal=`calcfits_bg test_image_real.fits = $BLINK_TEST_DATADIR/eda2/test_image_real.fits  | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi


echo "calcfits_bg test_image_imag.fits = $BLINK_TEST_DATADIR/eda2/test_image_imag.fits"
calcfits_bg test_image_imag.fits = $BLINK_TEST_DATADIR/eda2/test_image_imag.fits
equal=`calcfits_bg test_image_imag.fits = $BLINK_TEST_DATADIR/eda2/test_image_imag.fits  | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

echo "Exiting script test_eda2_data_20220118_ch204.sh with exit code = $exit_code"
exit $exit_code
