#!/bin/bash

#SBATCH -p EM
#SBATCH -J io-test
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --time=5:00:00 

DATADIR=/ocean/projects/cis240100p/gwilkins

datasets=(
    s3d/stat_planar.1.1000E-03.field.d64
    nyx/temperature.f32
    hacc/vx.f32
    cesm/V_1_26_1800_3600.f32
)

compressors=(
    sz
)

error_bounds=(
    0.1
    0.01
    0.001
    0.0001
    0.00001
    0.000001
)

io_methods=(
    netcdf
    hdf5
)

cd /jet/home/gwilkins/compression-energy/src/io_profiling
make
for d in ${datasets[@]}; do
cp $DATADIR/$d $LOCAL
for c in ${compressors[@]}; do
for eb in ${error_bounds[@]}; do
for io in ${io_methods[@]}; do
    ./io_test $c $d $eb $io
done
done
done
rm $LOCAL/$d
done
