#!/bin/bash

#SBATCH -J parallel-io          # Job name
#SBATCH -o parallel-io.o%j       # Name of stdout output file
#SBATCH -e parallel-io.e%j       # Name of stderr error file
#SBATCH -p skx             # Queue (partition) name
#SBATCH -N 24               # Total # of nodes 
#SBATCH -n 1024             # Total # of mpi tasks
#SBATCH -t 1:00:00        # Run time (hh:mm:ss)

datasets=(
   nyx/temperature.f32
)

compressors=(
    sz
    sz3
    zfp
    qoz
)

error_bounds=(
    0.001
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
