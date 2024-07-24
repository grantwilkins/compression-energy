#!/bin/bash

#SBATCH -J serial-compress

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=6:00:00 
#SBATCH --gres=gpu:1



module load amd-uprof
cd /home/gwilkins/compression-energy/src/compress_decompress
mkdir -p ./serial-compress/
make
AMDuProfCLI timechart --event power --interval 100 --duration 99999 -o ./serial-compress ./compress_cpu

