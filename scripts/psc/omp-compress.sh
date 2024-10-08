#!/bin/bash

#SBATCH -p EM
#SBATCH -J omp-compress
#SBATCH --nodes=1
#SBATCH -n 72
#SBATCH --time=5:00:00


datasets=(
        cesm/V_1_26_1800_3600.f32
 )

compressors=(
    sz3
    zfp
)

error_bounds=(
    0.001
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
