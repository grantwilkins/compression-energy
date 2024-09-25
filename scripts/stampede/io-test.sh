#!/bin/bash

#SBATCH -J io-test          # Job name
#SBATCH -o io-test.o%j       # Name of stdout output file
#SBATCH -e io-test.e%j       # Name of stderr error file
#SBATCH -p spr             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1             # Total # of mpi tasks
#SBATCH -t 6:00:00        # Run time (hh:mm:ss)

datasets=(
    nyx/temperature.f32
    s3d/stat_planar.1.1000E-03.field.d64
    cesm/V_1_26_1800_3600.f32
    hacc/vx.f32
)

compressors=(
    None
)

error_bounds=(
    0.0
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
