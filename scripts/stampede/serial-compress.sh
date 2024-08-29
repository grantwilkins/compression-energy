#!/bin/bash

#SBATCH -J serial-compress           # Job name
#SBATCH -o serial-compress.o%j       # Name of stdout output file
#SBATCH -e serial-compress.e%j       # Name of stderr error file
#SBATCH -p spr             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1             # Total # of mpi tasks
#SBATCH -t 6:00:00        # Run time (hh:mm:ss)

datasets=(
    nyx/temperature.f32
    stat_planar.1.1000E-03.field.d64
    cesm/V_1_26_1800_3600.f32
    hacc/vx.f32
)

compressors=(
    sz
    sz3
    zfp
    qoz
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
for j in ${compressors[@]}; do
for k in ${error_bounds[@]}; do
    ./compress_cpu $j $i $k
done
done
done
