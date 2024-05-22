#!/bin/bash

# /home/msok/Desktop/PAWSEY/PaCER/logbook/20230908_testing_BLINK_PIPELINE_cable_and_geom_correction_applycal.odt
# see /media/msok/80f59a5b-8bba-4392-8b31-a4f02fbee2f1/mwa/vcs/1276619416/1276619416_new/cotter/cotter_only_geomandcablecorr_applycal/BLINK_PIPELINE

n_channels=768
start_coarse_channel=133
fine_bw=0.04

module_path=`which module`

if [[ $PAWSEY_CLUSTER = "setonix" ]]; then
   pipeline_module=blink-pipeline/maingpu
   
   module reset 
   module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0 /software/setonix/2023.08/modules/zen3/gcc/12.2.0/libraries
   module load msfitslib/devel   
   module use /software/projects/director2183/msok/setonix/2023.08/modules/zen3/gcc/12.2.0/ /software/projects/director2183/setonix/2023.08/modules/zen3/gcc/12.2.0
   module load ${pipeline_module} blink_test_data/devel
fi
# salloc --partition gpuq --time 1:00:00 --nodes=1

echo "ln -sf $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch133.dat"
ln -sf $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch133.dat

echo "ln -sf $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch156.dat"
ln -sf $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch156.dat

if [[ ! -s 20200619163000.metafits ]]; then
   echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits ."
   cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits .
else
   echo "INFO : file 20200619163000.metafits already exists"
fi   

ch=0
while [[ $ch -lt 32 ]];
do
   ch_str=`echo $ch | awk '{printf("%03d",$1);}'`
   echo "cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan${ch_str}_xx.txt ."
   cp $BLINK_TEST_DATADIR/mwa/1276619416/calsolutions_chan${ch_str}_xx.txt .
   
   mkdir -p ch${ch_str}
   
   ch=$(($ch+1))
done   

echo "------------------- preparation completed -------------------"

ch=0
freq_mhz=`echo $ch | awk -v start_coarse_channel=${start_coarse_channel} -v fine_bw=${fine_bw} '{freq_mhz=start_coarse_channel*1.28;printf("%.4f\n",freq_mhz);}'`
coarse_channel=$start_coarse_channel

echo "blink_pipeline -c 4 -C -1 -t 1.00s -o ch -n 8192 -f ${freq_mhz} -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 100 -r -L -G -s calsolutions -r -V 100 -A 21,25,58,71,80,81,92,101,108,114,119,125 1276619416_1276619418_ch${coarse_channel}.dat > 133.out 2>&1"
blink_pipeline -c 4 -C -1 -t 1.00s -o ch -n 8192 -f ${freq_mhz} -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 100 -r -L -G -s calsolutions -r -V 100 -A 21,25,58,71,80,81,92,101,108,114,119,125 1276619416_1276619418_ch${coarse_channel}.dat > 133.out 2>&1

echo "------------------- blink_pipeline execution completed -------------------"

# Running comparsison using existing script (gzipped FITS files - hence 1):
# WARNING : needs a bit of work to use compare_blink_images_gpu_vs_cpu.sh to compare GPU version with the template files produced by the CPU version 
#           as this is slightly different scenario and FITS files are not in one-to-one relation
echo "../experiments/compare_one_to_one.sh ${BLINK_TEST_DATADIR}/mwa/1276619416/results/full_imaging/ "" ch 1 0 > comparison.out 2>&1"
../experiments/compare_one_to_one.sh ${BLINK_TEST_DATADIR}/mwa/1276619416/results/full_imaging/ "" ch 1 0 > comparison.out 2>&1 
equal_count=`grep COMPARISON comparison.out  | grep EQUAL | wc -l`
if [[ $equal_count -eq 36 ]]; then
   echo "OK result of the test"
   exit_code=0
else
   echo "WARNING : test failed number of EQUAL files = $equal_count vs. expected 54 !"
   exit_code=1
fi

me=$(basename "$0")
echo "Exiting script $me with exit code = $exit_code"
exit $exit_code

echo "------------------- full test execution completed -------------------"

