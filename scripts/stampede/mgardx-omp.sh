#!/bin/bash

#SBATCH -J mgardx-omp           # Job name
#SBATCH -o mgardx-omp.o%j       # Name of stdout output file
#SBATCH -e mgardx-omp.e%j       # Name of stderr error file
#SBATCH -p icx             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1             # Total # of mpi tasks
#SBATCH -t 6:00:00        # Run time (hh:mm:ss)

datasets=(
    nyx/temperature.f32
    cesm/V_1_26_1800_3600.f32
)

error_bounds=(
    0.001
)

cd $HOME/compression-energy/src/intel_compress_decompress/
make
CORES=(1 2 4 8 16 32 64)
for i in ${datasets[@]}; do
for k in ${error_bounds[@]}; do
for c in ${CORES[@]}; do
    export OMP_NUM_THREADS=$c
    ./build/mgardx_omp $i $k
done
done
done
