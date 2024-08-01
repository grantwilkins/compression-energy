# Compression Energy Experiments

This repository is the implementation of a testbed to characterize the energy efficiency of scientific, floating-point lossy compressors. In particular, we are interested in demonstrating that lossy compression can be used to not just reduce data storage sizes, but also reduce the energy required to store and transmit data.

The following are the main components of this repository:

- `src/compress_decompress`: This directory contains the implementation of the compress and decompress experiments. In this experiment, we will compress and decompress several datasets using different lossy compression algorithms and compare the results.
- `src/hdf5_transfer`: This directory contains the implementation of transmitting data to a parallel file system using HDF5. We will use this experiment to measure the energy required to transfer data to a parallel file system.
- `scripts`: This directory contains the scripts to run the experiments and collect the power draw across all CPU cores.
- `analysis`: This directory contains the scripts to generate the plots from the experiment results.

## Required Software Packages

This repository requires the following software packages:

- [Libpressio](https://github.com/robertu94/libpressio) to implement several lossy compressors.
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) to transfer the compressed data.
- [AMDuProf](https://developer.amd.com/amd-uprof/) to collect the power draw across all CPU cores.
- [NVIDIA-smi](https://developer.nvidia.com/nvidia-system-management-interface) to collect the power draw of the GPU.


## Variables

- Datasets:
  - Miranda: 3072x3072x3072, 106GB, Single Precision, 1 field
  - HACC: 1073726487, 19GB, Single Precision, 6 fields
  - S3D: 500x500x500, 44GB, Double Precision, 11 fields
  - NYX: 512x512x512, 4GB, Single Precision, 6 fields
- Single-Threaded Compressors:
  - SZ2
  - SZ3
  - ZFP
  - SZx
  - MGARD
  - MGARDx
  - ZFP-X
- Multi-Threaded Compressors:
  - SZ-OMP 
  - SZ3-OMP
  - ZFP-OMP
  - SZx-OMP
  - MGARD-OMP
- GPU Compressors:
  - FZ-GPU
  - SZ-GPU
  - ZFP-GPU
  - SZx-GPU
  - MGARD-GPU
  - cuSZp
