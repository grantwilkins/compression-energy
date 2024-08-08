#!/bin/bash

#SBATCH -J serial-compress
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=12:00:00 
#SBATCH --gres=gpu:1

datasets=(
    hacc/vz.f32
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

module load amd-uprof/4.2.850
module unload amd-uprof/4.2.850
module load amd-uprof/4.2.850
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
