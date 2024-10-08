#!/bin/bash

#SBATCH -p EM
#SBATCH -J serial-compress
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --time=6:00:00 

datasets=(
    hacc/vy.f32
)

compressors=(
    qoz
)

error_bounds=(
    0.1
    0.01
    0.001
    0.0001
    0.00001
    0.000001
)

cd /jet/home/gwilkins/compression-energy/src/intel_compress_decompress
make
for i in ${datasets[@]}; do
for j in ${compressors[@]}; do
for k in ${error_bounds[@]}; do
    ./compress_cpu $j $i $k
done
done
done
