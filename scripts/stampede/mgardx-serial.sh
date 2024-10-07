#!/bin/bash

#SBATCH -J mgard-serial           # Job name
#SBATCH -o mgard-serial.o%j       # Name of stdout output file
#SBATCH -e mgard-serial.e%j       # Name of stderr error file
#SBATCH -p icx             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1             # Total # of mpi tasks
#SBATCH -t 6:00:00        # Run time (hh:mm:ss)

datasets=(
    nyx/baryon_density.f32
    cesm/V_1_26_1800_3600.f32
)

error_bounds=(
    0.1
    0.01
    0.001
    0.0001
    0.00001
    0.000001
)

cd $HOME/compression-energy/src/intel_compress_decompress/
make
for i in ${datasets[@]}; do
for k in ${error_bounds[@]}; do
    ./build/mgardx $i $k
done
done

