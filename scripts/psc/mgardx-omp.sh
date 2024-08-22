#!/bin/bash

#SBATCH -p EM
#SBATCH -J mgardx-omp
#SBATCH --nodes=1
#SBATCH -n 72
#SBATCH --time=6:00:00 

datasets=(
    hacc/vx.f32
)


error_bounds=(
    0.001
)

CORES=(32 64)
cd /jet/home/gwilkins/compression-energy/src/intel_compress_decompress
for d in ${datasets[@]}; do
for eb in ${error_bounds[@]}; do
for c in ${CORES[@]}; do
    export OMP_NUM_THREADS=$c
    ./build/mgardx_omp $d $eb
done
done
done
