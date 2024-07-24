#!/bin/bash

#SBATCH -J omp-compress

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=6:00:00 
#SBATCH --gres=gpu:1


CORES = (1 2 4 8 16 32 64 128)
module load amd-uprof
cd /home/gwilkins/compression-energy/src/compress_decompress
mkdir -p ./omp-compress/
make
for i in ${CORES[@]}; do
    AMDuProfCLI timechart --event power --interval 100 --duration 99999 -o ./omp-compress ./compress_cpu_omp $i
done

