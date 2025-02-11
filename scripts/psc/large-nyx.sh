#!/bin/bash

#SBATCH -p EM
#SBATCH -J large-nyx
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --time=6:00:00 

datasets=(
    nyx/temperature_2.f32
    nyx/temperature_3.f32
    nyx/temperature_4.f32
    nyx/temperature_5.f32
)

compressors=(
    sz
    sz3
    qoz
    zfp
)

error_bounds=(
    0.1
    0.01
    0.001
    0.0001
    0.00001
)

cd /jet/home/gwilkins/compression-energy/src/intel_compress_decompress
make
for i in ${datasets[@]}; do
for j in ${compressors[@]}; do
for k in ${error_bounds[@]}; do
    ./compress_cpu_large_nyx $j $i $k
done
done
done
