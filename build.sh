#!/bin/bash -e

build_type="cpu"
build_dir=build
if [[ -n "$1" && "$1" != "-" ]]; then
   build_type=$1
fi

echo "PARAMETERS:"
echo "build_dir = $build_dir"

cmake_options=""
# DEBUG : -DCMAKE_BUILD_TYPE=Debug
if [[ $build_type == "gpu" ]]; then
   build_dir=build_gpu
fi
if [[ $build_type == "gpu" ]]; then
   # WARNING : -DCUDA_SM=70 is specific for Topaz, the default is sm=61 for msok's laptop, but laptop/desktop is compiled using normal cmake 
   # cmake_options="-DUSE_HIP=ON -DCUDA_SM=70"
   cmake_options="-DUSE_HIP=ON -DCMAKE_BUILD_TYPE=Debug"
   echo "INFO (build.sh) : build type ${build_type} detected -> added cmake_options = ${cmake_options}"
else
   echo "INFO (build.sh) : CPU version"
fi


# First, you need to source the bash library
module load bash-utils
source "${BASH_UTILS_DIR}/build_utils.sh"


PROGRAM_NAME=blink-pipeline
PROGRAM_VERSION=main

 
# the following function sets up the installation path according to the
# cluster the script is running on and the first argument given. The argument
# can be:
# - "group": install the software in the group wide directory
# - "user": install the software only for the current user
# - "test": install the software in the current working directory 
process_build_script_input user #group 


# load all the modules required for the program to compile and run.
# the following command also adds those module names in the modulefile
# that this script will generate.
echo "Loading required modules ..."
if [ $PAWSEY_CLUSTER = "setonix" ]; then
   module reset
   module load cmake/3.24.3
   if [[ $build_type == "gpu" ]]; then
      print_run module_load blink_test_data/devel  blink_astroio/master  blink_correlator/master blink-imager-gpu/devel blink_preprocessing/main rocm/5.4.3
   else
      print_run module_load blink_test_data/devel  blink_astroio/master  blink_correlator/master blink-imager-cpu/devel blink_preprocessing/main rocm/5.4.3
   fi
else
   print_run module purge
   print_run module_load gcc/8.3.0 cascadelake blink_test_data/devel blink-correlator/devel-cmplx blink-imager/cristian-dev
fi
# cmake is only required at build time, so we use the normal module load
# build your software..
echo "Building the software.."

[ -d ${build_dir} ] || mkdir ${build_dir}
cd ${build_dir}
# Turns out we need to compile with HIPCC if AstroIO and the other libraries were compiled with GPU support. This is because
# Voltages and Visibilities classes derive from MemoryBuffer, a template class in a header file make use of GPU calls.
print_run cmake .. -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_BUILD_TYPE=Debug
make VERBOSE=1
# make test
# Install the software
# make install

# echo "Create the modulefile.."
# create_modulefile

echo "Done."


