CC=gcc
MPICC=mpicc

HDF5_INC=$(TACC_HDF5_INC)
HDF5_LIB=$(TACC_HDF5_LIB)
LDLIBS=$(shell pkg-config --libs libpressio papi netcdf pnetcdf) -D_POSIX_C_SOURCE=199309L
CFLAGS=-std=c99 -lm -g -lpapi -fopenmp $(shell pkg-config --cflags libpressio papi netcdf pnetcdf)

ALL=io_test
.PHONY: all clean

io_test: io_test.c
	$(MPICC) $(CFLAGS) -o io_test io_test.c -I$(HDF5_INC) -L$(HDF5_LIB) $(LDLIBS) -lhdf5 -lmpi

clean:
	rm io_test
