#!/bin/bash

#!/bin/bash -e 
#SBATCH --account=mwavcc
#SBATCH --partition=mwa
#SBATCH --time=23:59:59
#SBATCH --export=NONE
#SBATCH --output=slurm-%A.out 
#SBATCH --error=./slurm-%j.err
# A nice way to execute a command and print the command line
function print_run {
        echo "Executing $@"
        echo "$@" > doit.sh 
        $@
}
# module use /software/setonix/unsupported/
# module use /software/projects/pawsey1154/cdipietrantonio/setonix/2025.08/modules/zen3/gcc/14.2.0/blink-pipeline-gpu/
module load blink-pipeline-gpu/main


export LD_LIBRARY_PATH=/software/projects/pawsey1154/msok/github/branches/msok-totalpower/pipeline/build/:$LD_LIBRARY_PATH
exec_dir=/software/projects/pawsey1154/msok/github/branches/msok-totalpower/pipeline/build/
if [[ -n "$1" && "$1" != "-" ]]; then
   exec_dir="$1"
fi

exec=${exec_dir}/test_totalpower

ls dyna*.fits > fits_list

print_run srun $exec fits_list 

for total_power_file in `ls *.total_power`
do
   exclude_ranges_file=${total_power_file%%total_power}exclude_ranges
   print_run srun ${exec_dir}/exclude_ranges_totalpower $total_power_file -o $exclude_ranges_file
done
   
