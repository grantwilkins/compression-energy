#!/bin/bash

#SBATCH -J serial-compress
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=12:00:00 
#SBATCH --gres=gpu:1

datasets=(
    s3d/stat_planar.1.1000E-03.field.d64
    nyx/temperature.f32
    hacc/vx.f32
    miranda/density.f32
)

compressors=(
    sz
    sz3
    zfp
    mgard
)

error_bounds=(
    0.1
    0.01
    0.001
    0.0001
    0.00001
    0.000001
)

module load amd-uprof/4.1.424
module unload amd-uprof/4.1.424
module load amd-uprof/4.1.424
cd /home/ac.gwilkins/compression-energy/src/compress_decompress
make
for i in ${datasets[@]}; do
for j in ${compressors[@]}; do
for k in ${error_bounds[@]}; do
    mkdir -p ./serial-compress/$i-$j-$k
    AMDuProfCLI timechart --event power --interval 500 --duration 99999 -o ./serial-compress/$i-$j-$k ./compress_cpu $j $i $k
done
done
done
