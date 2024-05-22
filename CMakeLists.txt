cmake_minimum_required(VERSION 3.10.0)
project(pipeline)
set(CMAKE_CXX_STANDARD 17)

include(CTest)

if(USE_HIP)
    message("INFO : HIP/GPU version")
    add_definitions("-DIMAGER_HIP")
else()
    message("INFO : CPU version")
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

# Pipeline library
file(GLOB pipeline_sources "src/*.cpp")
file(GLOB pipeline_headers "src/*.hpp")

add_library(pipeline SHARED ${pipeline_sources})
set_target_properties(pipeline PROPERTIES PUBLIC_HEADER "${pipeline_headers}")
target_link_libraries(pipeline ${CORRELATIONLIB} ${BLINK_ASTROIO} ${BLINK_PREPROCESSING} ${PACERIMAGERLIB} ${MSFITSLIB} ) #${LIBPAL_LIB})

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(pipeline OpenMP::OpenMP_CXX)
endif()


install(TARGETS pipeline
   LIBRARY DESTINATION "lib"
   PUBLIC_HEADER DESTINATION "include"
)

add_definitions(-D_UNIX -D_PACER_PROFILER_ON_)


# APPS
add_executable(blink_pipeline apps/blink_pipeline.cpp)
target_link_libraries(blink_pipeline pipeline stdc++fs)
install(TARGETS blink_pipeline DESTINATION "bin")

# TESTS

# Proper test using the same program (no special test-program)
find_program (BASH_PROGRAM bash)
# add_test(NAME test_eda2_data_20220118_ch204.sh COMMAND ../tests/test_eda2_data_20220118_ch204.sh)
# add_test(NAME test_mwa_obsid1276619416_full_imaging.sh COMMAND ../tests/test_mwa_obsid1276619416_full_imaging.sh)
