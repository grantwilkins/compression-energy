#!/bin/bash

#SBATCH -J serial-compress
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=6:00:00 
#SBATCH --gres=gpu:1

datasets=(
    nyx/temperature.f32
    nyx/density.f32
    nyx/velocity_x.f32
    nyx/velocity_y.f32
    nyx/velocity_z.f32
    hacc/vx.f32
    hacc/vy.f32
    hacc/vz.f32
    hacc/xx.f32
    hacc/yy.f32
    hacc/zz.f32
    miranda/density.f32
)

compressors=(
    sz
    sz3
    zfp
    mgard
)

module load amd-uprof/4.1.424
module unload amd-uprof/4.1.424
module load amd-uprof/4.1.424
cd /home/ac.gwilkins/compression-energy/src/compress_decompress
make
for i in ${datasets[@]}; do
for j in ${compressors[@]}; do
    mkdir -p ./serial-compress/$i-$j
    AMDuProfCLI timechart --event power --interval 100 --duration 99999 -o ./serial-compress/$i-$j ./compress_cpu $j $i
done
done
