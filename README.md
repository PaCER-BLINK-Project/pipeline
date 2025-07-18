# BLINK Pipeline

This repository implements the BLINK imaging pipeline for FRB searches.

## Dependencies

You will need the following libraries to compile the project:

- `astroio`, `correlation`, `preprocessing`, `imager` part of the BLINK project.

Other dependencies are optional:

- OpenMP
- HIP/ROCm for AMD GPU acceleration.
- CUDA for NVIDIA GPU acceleration.

## Compiling

The compilation process is handled by CMake. To build a CPU only version, do the following:

```
mkdir build; cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make 
```

A CUDA enabled build can be obtained by simply adding the `-DUSE_CUDA=ON` argument to `cmake`.

To compile the code with HIP support, you will need to specify `-DUSE_HIP=ON -DCMAKE_CXX_COMPILER=hipcc`. 

To run tests, execute `make test`.

Available CMake flags are:

- `USE_OPENMP` (default: `ON`): run the correlation function in parallel on multiple cores.
- `USE_HIP` (default: `OFF`): build the library using HIP to enable AMD GPU support.
- `USE_CUDA` (default: `OFF`): build the library using CUDA to enable NVIDIA GPU support.