#!/bin/bash

#SBATCH -J io-test          # Job name
#SBATCH -o io-test.o%j       # Name of stdout output file
#SBATCH -e io-test.e%j       # Name of stderr error file
#SBATCH -p skx             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1             # Total # of mpi tasks
#SBATCH -t 6:00:00        # Run time (hh:mm:ss)

datasets=(
    s3d/stat_planar.1.1000E-03.field.d64
    hacc/vx.f32
)

compressors=(
    None
)

error_bounds=(
    0.1
    0.01
    0.001
    0.0001
    0.00001
    0.000001
)

cd $HOME/compression-energy/src/io_profiling
make
CORES=(1 2 4 8 16 32 64)
for d in ${datasets[@]}; do
for eb in ${error_bounds[@]}; do
for c in ${compressors[@]}; do
	./io_test $c $d $eb
done
done
done
