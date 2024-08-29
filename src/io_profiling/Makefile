CC=gcc
MPICC=mpicc
LDLIBS=$(shell pkg-config --libs libpressio papi)
CFLAGS=-std=c99 -lm -g -lpapi $(shell pkg-config --cflags libpressio papi)	

IO_INCLUDE=-I${TACC_HDF5_INC} -I${TACC_NETCDF_INC} -I${TACC_PNETCDF_INC}
IO_CFLAGS=-std=c99 -lm -fopenmp -g -lpapi -lhdf5 -lnetcdf -lpnetcdf $(shell pkg-config --cflags libpressio papi)
ALL=io_test
.PHONY: all clean

prepare_compressed_files: prepare_compressed_files.c
	$(CC) $(CFLAGS) -o prepare_compressed_files prepare_compressed_files.c $(LDLIBS)

io_test: io_test.c
	$(CC) $(CFLAGS) -o io_test io_test.c $(LDLIBS)


