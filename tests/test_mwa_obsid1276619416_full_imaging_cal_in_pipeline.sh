#!/bin/bash -e

# No module loads here! Modules are loaded in the build script or from the terminal, or other env scripts.

# A nice way to execute a command and print the command line
function print_run {
	echo "Executing $@"
	$@
}

# Create a test directory, so test files do no pollute the build directory. 
[ -e test_full_imaging ] || mkdir test_full_imaging
cd test_full_imaging

print_run ln -sf $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch133.dat

print_run ln -sf $BLINK_TEST_DATADIR/mwa/1276619416/voltages/1276619416_1276619418_ch156.dat

if [[ ! -s 20200619163000.metafits ]]; then
   print_run cp $BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits .
else
   echo "INFO : file 20200619163000.metafits already exists"
fi   

SOL_FILE=${BLINK_TEST_DATADIR}/mwa/1276619416/1276625432.bin

print_run mkdir -p ch000 ch001 ch750

print_run ../blink_pipeline -c 4 -C 0 -t 1.00s -o ch000/ -n 8192 -f 169.6000 -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 0 -r -L -G -s ${SOL_FILE} -b 0 -r -V 100 -A 21,25,58,71,80,81,92,101,108,114,119,125   1276619416_1276619418_ch133.dat 
print_run ../blink_pipeline -c 4 -C 1 -t 1.00s -o ch001/ -n 8192 -f 169.6400 -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 0 -r -L -G -s ${SOL_FILE} -b 0 -r -V 100 -A 21,25,58,71,80,81,92,101,108,114,119,125  1276619416_1276619418_ch133.dat 
print_run ../blink_pipeline -c 4 -C 14 -t 1.00s -o ch750/ -n 8192 -f 199.6000 -F 30 -M 20200619163000.metafits -U 1592584240 -w N -v 0 -r -L -G -s ${SOL_FILE} -b 23 -r -V 100 -A 21,25,58,71,80,81,92,101,108,114,119,125  1276619416_1276619418_ch156.dat 


echo "Checking correctness of results.."


for CH_STR in ch000 ch001 ch750; 
do
cd ${CH_STR}
for fitsfile in `ls -1 *.fits`;
do
calcfits_bg $fitsfile = $BLINK_TEST_DATADIR/mwa/1276619416/results/full_imaging/${CH_STR}/$fitsfile.gz
done
cd ..
done

