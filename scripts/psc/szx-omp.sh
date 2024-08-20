#!/bin/bash

#SBATCH -p EM
#SBATCH -J szx-omp
#SBATCH --nodes=1
#SBATCH -n 72
#SBATCH --time=6:00:00 

datasets=(
    nyx/temperature.f32
    cesm/V_1_26_1800_3600.f32
    s3d/stat_planar.1.1000E-03.field.d64
    hacc/vx.f32
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
for eb in ${error_bounds[@]}; do
for c in ${CORES[@]}; do
    export OMP_NUM_THREADS=$c
    ./szx_omp $d $eb
done
done
done
