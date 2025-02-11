#!/bin/bash

#SBATCH -p EM
#SBATCH -J large-nyx
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --time=6:00:00 


cd /jet/home/gwilkins/compression-energy/src/intel_compress_decompress
make
./inflate_nyx /ocean/projects/cis240100p/gwilkins/nyx/temperature.f32 64 /ocean/projects/cis240100p/gwilkins/nyx/temperature_4.f32
