#!/bin/bash

#SBATCH -J serial-compress
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=12:00:00 
#SBATCH --gres=gpu:1

datasets=(
    hacc/vx.f32
    hacc/xx.f32
    miranda/density.f32
)

compressors=(
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
