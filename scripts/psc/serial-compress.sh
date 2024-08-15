#!/bin/bash

#SBATCH -p EM
#SBATCH -J serial-compress
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --time=3:00:00 

datasets=(
    s3d/stat_planar.1.1000E-03.field.d64
)

compressors=(
    mgard
)

error_bounds=(
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
