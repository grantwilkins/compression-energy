#!/bin/bash

#SBATCH -p EM
#SBATCH -J mgardx-compress
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --time=12:00:00 

datasets=(
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

cd /jet/home/gwilkins/compression-energy/src/intel_compress_decompress
for d in ${datasets[@]}; do
for eb in ${error_bounds[@]}; do
    ./build/mgardx $d $eb
done
done
