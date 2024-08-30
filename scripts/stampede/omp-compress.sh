#!/bin/bash

#SBATCH -J omp-compress           # Job name
#SBATCH -o omp-compress.o%j       # Name of stdout output file
#SBATCH -e omp-compress.e%j       # Name of stderr error file
#SBATCH -p spr             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1             # Total # of mpi tasks
#SBATCH -t 12:00:00        # Run time (hh:mm:ss)

datasets=(
    s3d/stat_planar.2.3500E-03.field.d64
)

compressors=(
    sz_omp
    sz3
    zfp
    qoz
)

error_bounds=(
    0.001
)

CORES=(1 2 4 8 16 32 64)
cd $HOME/compression-energy/src/intel_compress_decompress/
make
for i in ${datasets[@]}; do
for j in ${compressors[@]}; do
for k in ${error_bounds[@]}; do
for c in ${CORES[@]}; do
     export OMP_NUM_THREADS=$c
    ./compress_cpu_omp $j $i $k
done
done
done
done
