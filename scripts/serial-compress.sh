#!/bin/bash

#SBATCH -J serial-compress

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=6:00:00 
#SBATCH --gres=gpu:1

datasets=(
    s3d/stat_planar.1.1000E-03.field.d64
    s3d/stat_planar.1.7000E-03.field.d64
    s3d/stat_planar.2.3500E-03.field.d64
    s3d/stat_planar.2.9000E-03.field.d64
    s3d/stat_planar.2.9950E-03.field.d64
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

module load amd-uprof
cd /home/ac.gwilkins/compression-energy/src/compress_decompress
make
for i in ${datasets[@]}; do
for j in ${compressors[@]}; do
    mkdir -p ./serial-compress/$i-$j
    AMDuProfCLI timechart --event power --interval 100 --duration 99999 -o ./serial-compress/$i-$J ./compress_cpu $j $i
done
done