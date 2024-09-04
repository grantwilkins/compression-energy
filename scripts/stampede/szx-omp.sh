#!/bin/bash

#SBATCH -J szx-omp         # Job name
#SBATCH -o szx-omp.o%j       # Name of stdout output file
#SBATCH -e szx-omp.e%j       # Name of stderr error file
#SBATCH -p skx           # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1             # Total # of mpi tasks
#SBATCH -t 6:00:00        # Run time (hh:mm:ss)

datasets=(
    hacc/vx.f32
    s3d/stat_planar.1.1000E-03.field.d64
    nyx/temperature.f32
    cesm/U_1_26_1800_3600.f32
)

error_bounds=(
    0.001
)

CORES=(1 2 4 8 16 32)
cd $HOME/compression-energy/src/intel_compress_decompress/
make
for i in ${datasets[@]}; do
for k in ${error_bounds[@]}; do
for c in ${CORES[@]}; do
    export OMP_NUM_THREADS=$c
    ./szx_omp $i $k
done
done
done
