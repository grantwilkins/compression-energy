#!/bin/bash

#SBATCH -p EM
#SBATCH -J szx-compress
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --time=6:00:00 

datasets=(
    nyx/temperature_2.f32
    nyx/temperature_3.f32
    nyx/temperature_4.f32
    nyx/temperature_5.f32
)

error_bounds=(
    0.001
)

cd /jet/home/gwilkins/compression-energy/src/intel_compress_decompress
make
for d in ${datasets[@]}; do
for eb in ${error_bounds[@]}; do
    ./szx_serial $d $eb
done
done
