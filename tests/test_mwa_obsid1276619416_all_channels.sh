#!/bin/bash -e

# TODO: this is just a base template script - currently output is nonsense but it shows the
# pipeline is able to handle an entire second of data. We need to fix the last part of the
# pipeline to allow the imaging software to process multiple coarse channels.

# No module loads here! Modules are loaded in the build script or from the terminal, or other env scripts.

# A nice way to execute a command and print the command line
function print_run {
	echo "Executing $@"
	$@
}

# Create a test directory, so test files do no pollute the build directory. 
[ -e test_all_channels ] || mkdir test_all_channels
cd test_all_channels

INPUT_FILES="/scratch/mwavcs/msok/1276619416/combined/1276619416_1276619418_ch1*.dat"
METAFITS="$BLINK_TEST_DATADIR/mwa/1276619416/20200619163000.metafits"
SOL_FILE=${BLINK_TEST_DATADIR}/mwa/1276619416/1276625432.bin


../blink_pipeline -c 4 -C -1 -t 1.00s -o ${MYSCRATCH}/1276619416_1276619418_images -n 8192 -f -1 -F 30 -M ${METAFITS} -U 1592584240 -w N -v 0 -r -L -G -s ${SOL_FILE} -b 0 -r -V 100 -A 21,25,58,71,80,81,92,101,108,114,119,125 ${INPUT_FILES} 




