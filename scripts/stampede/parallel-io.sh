#!/bin/bash

#SBATCH -J parallel-io          # Job name
#SBATCH -o parallel-io.o%j       # Name of stdout output file
#SBATCH -e parallel-iot.e%j       # Name of stderr error file
#SBATCH -p spr             # Queue (partition) name
#SBATCH -N 2               # Total # of nodes 
#SBATCH -n 150             # Total # of mpi tasks
#SBATCH -t 6:00:00        # Run time (hh:mm:ss)

datasets=(
   nyx/temperature.f32
)

compressors=(
    zfp
)

error_bounds=(
    0.1
    0.01
    0.001
    0.0001
    0.00001
    0.000001
)

cd $HOME/compression-energy/src/parallel_io
make
CORES=(1 2 4 8 16 32 64)
for d in ${datasets[@]}; do
for eb in ${error_bounds[@]}; do
for c in ${compressors[@]}; do
	ibrun ./parallel_io $c $d $eb
done
done
done
