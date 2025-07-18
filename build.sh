#!/bin/bash -e


# First, you need to source the bash library
module load bash-utils
source "${BASH_UTILS_DIR}/build_utils.sh"


PROGRAM_NAME=blink-pipeline-gpu
PROGRAM_VERSION=mvded

 
# the following function sets up the installation path according to the
# cluster the script is running on and the first argument given. The argument
# can be:
# - "group": install the software in the group wide directory
# - "user": install the software only for the current user
# - "test": install the software in the current working directory 
process_build_script_input user


# load all the modules required for the program to compile and run.
# the following command also adds those module names in the modulefile
# that this script will generate.
echo "Loading required modules ..."
module reset
module load cmake/3.30.5
print_run module_load blink_test_data/devel blink-astroio/master blink-correlator/master blink-imager-gpu/marcin-branch blink-preprocessing/main blink-dedispersion/main rocm/6.4.1


# cmake is only required at build time, so we use the normal module load
# build your software..
echo "Building the software.."
build_dir=build
[ -d ${build_dir} ] || mkdir ${build_dir}
cd ${build_dir}
# Turns out we need to compile with HIPCC if AstroIO and the other libraries were compiled with GPU support. This is because
# Voltages and Visibilities classes derive from MemoryBuffer, a template class in a header file make use of GPU calls.
print_run cmake .. -DUSE_OPENMP=ON -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DUSE_HIP=ON -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_C_COMPILER=hipcc -DCMAKE_CXX_FLAGS=-O0 -DCMAKE_BUILD_TYPE=Debug
make VERBOSE=1


# Install the software
make install

echo "Create the modulefile.."
create_modulefile

echo "Done."


