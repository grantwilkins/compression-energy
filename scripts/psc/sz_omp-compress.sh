#!/bin/bash

#SBATCH -J sz-omp-compress
#SBATCH -p EM
#SBATCH --nodes=1
#SBATCH -n 72
#SBATCH --time=12:00:00 


datasets=(
    s3d/stat_planar.1.3500E-03.field.d64
)

compressors=(
    sz_omp
)

error_bounds=(
    0.1
    0.01
    0.001
    0.0001
    0.00001
    0.000001
)

CORES=(1 2 4 8 16 32 64)
cd /jet/home/gwilkins/compression-energy/src/intel_compress_decompress
make
for d in ${datasets[@]}; do
for c in ${compressors[@]}; do
for e in ${error_bounds[@]}; do
for i in ${CORES[@]}; do
    export OMP_NUM_THREADS=$i
    ./compress_cpu_omp $c $d $e
done
done
done
done
