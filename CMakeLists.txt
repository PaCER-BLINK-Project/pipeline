cmake_minimum_required(VERSION 3.10.0)
project(blink_astroio)
set(CMAKE_CXX_STANDARD 17)
include(CTest)

option(USE_CUDA "Compile the code with NVIDIA GPU support." OFF)
option(USE_HIP "Compile the code with AMD GPU support." OFF)
option(USE_OPENMP "Compile with OpenMP enabled." ON)

if(USE_CUDA)
enable_language(CUDA CXX)
endif()

option(USE_OPENMP "Compile with OpenMP enabled." ON)

# External dependencies
find_library(CFITSIO_LIB cfitsio HINTS ENV LD_LIBRARY_PATH)
if( NOT CFITSIO_LIB )
    message(FATAL_ERROR "CFITSIO library not found.")
endif()

# find_library(LIBPAL_LIB pal PATHS ENV LD_LIBRARY_PATH)

find_library(MSFITSLIB NAMES libmsfitslib.so HINTS ENV LD_LIBRARY_PATH)
if(NOT MSFITSLIB)
   message(FATAL_ERROR "MSFITSLIB library not found")
endif()

find_library(PACERIMAGERLIB NAMES libimaging.so HINTS ENV LD_LIBRARY_PATH)
if(NOT PACERIMAGERLIB)
   message(FATAL_ERROR "PACERIMAGERLIB library not found")
endif()

find_library(CORRELATIONLIB correlation HINTS ENV LD_LIBRARY_PATH)
if (NOT CORRELATIONLIB)
    message(FATAL_ERROR "CORRELATION library not found.")
endif()

find_library(BLINK_ASTROIO blink_astroio HINTS ENV LD_LIBRARY_PATH)
if (NOT BLINK_ASTROIO)
   message(FATAL_ERROR "BLINK_ASTROIO library not found.")
endif()

find_library(BLINK_PREPROCESSING blink_preprocessing HINTS ENV LD_LIBRARY_PATH)
if (NOT BLINK_PREPROCESSING)
   message(FATAL_ERROR "BLINK_PREPROCESSING library not found.")
endif()

find_library(BLINK_DEDISP blink_dedispersion HINTS ENV LD_LIBRARY_PATH)
if (NOT BLINK_DEDISP)
   message(FATAL_ERROR "BLINK_DEDISP library not found.")
endif()

# Pipeline library
file(GLOB pipeline_sources "src/*.cpp")
file(GLOB pipeline_headers "src/*.hpp")
if(USE_CUDA OR USE_HIP)
   file(GLOB gpu_sources "src/gpu/*.cpp")
   file(GLOB gpu_headers "src/gpu/*.hpp")
   list(APPEND pipeline_sources ${gpu_sources})
endif()

if(USE_CUDA)
set_source_files_properties(${pipeline_sources} "apps/blink_pipeline.cpp" PROPERTIES LANGUAGE CUDA)
endif()

add_library(pipeline SHARED ${pipeline_sources})
set_target_properties(pipeline PROPERTIES PUBLIC_HEADER "${pipeline_headers}")
target_link_libraries(pipeline ${CORRELATIONLIB} ${BLINK_ASTROIO} ${BLINK_PREPROCESSING} ${PACERIMAGERLIB} ${MSFITSLIB} ${BLINK_DEDISP} )


install(TARGETS pipeline
   LIBRARY DESTINATION "lib"
   PUBLIC_HEADER DESTINATION "include"
)
install(FILES ${gpu_headers} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gpu")

# APPS
add_executable(blink_pipeline apps/blink_pipeline.cpp)
target_link_libraries(blink_pipeline pipeline stdc++fs)
install(TARGETS blink_pipeline DESTINATION "bin")

add_executable(test_totalpower apps/test_totalpower.cpp)
target_link_libraries(test_totalpower pipeline stdc++fs)
install(TARGETS test_totalpower DESTINATION "bin")

add_executable(exclude_ranges_totalpower apps/exclude_ranges_totalpower.cpp)
install(TARGETS exclude_ranges_totalpower DESTINATION "bin")


if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(pipeline OpenMP::OpenMP_CXX)
    target_link_libraries(blink_pipeline OpenMP::OpenMP_CXX)
endif()

# TESTS

# Proper test using the same program (no special test-program)
find_program (BASH_PROGRAM bash)
# add_test(NAME test_eda2_data_20220118_ch204.sh COMMAND ../tests/test_eda2_data_20220118_ch204.sh)
# add_test(NAME test_mwa_obsid1276619416_full_imaging.sh COMMAND ../tests/test_mwa_obsid1276619416_full_imaging.sh)
