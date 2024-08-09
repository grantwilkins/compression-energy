#!/bin/bash

#SBATCH -J omp-compress

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=12:00:00 
#SBATCH --gres=gpu:1


datasets=(
    hacc/vy.f32
    cesm/V_1_26_1800_3600.f32
)

compressors=(
    zfp
)

error_bounds=(
    0.001
)

CORES=(1 2 4 8 16 32 64)
module load amd-uprof/4.2.850
module unload amd-uprof/4.2.850
module load amd-uprof/4.2.850
cd /home/ac.gwilkins/compression-energy/src/compress_decompress
mkdir -p ./omp-compress/
make
for d in ${datasets[@]}; do
for c in ${compressors[@]}; do
for e in ${error_bounds[@]}; do
for i in ${CORES[@]}; do
    mkdir -p ./omp-compress/$d-$c-$e-$i
    export OMP_NUM_THREADS=$i
    AMDuProfCLI timechart --event power --interval 500 --duration 99999 -o ./omp-compress/$d-$c-$e-$i ./compress_cpu_omp $c $d $e
done
done
done
done
