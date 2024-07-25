#!/bin/bash

#SBATCH -J serial-compress

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=6:00:00 
#SBATCH --gres=gpu:1

datasets=(   
    hacc/vy.f32
    hacc/vz.f32
    hacc/xx.f32
    hacc/yy.f32
    hacc/zz.f32
    miranda/density.f32
)

compressors=(
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

