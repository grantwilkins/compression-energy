CC=gcc
MPICC=mpicc

SZX_DIR := /home1/10191/gfw/SZx/build
SZX_INC := $(SZX_DIR)/include
SZX_LIB := $(SZX_DIR)/lib

HDF5_INC=$(TACC_HDF5_INC)
HDF5_LIB=$(TACC_HDF5_LIB)
LDLIBS=$(shell pkg-config --libs libpressio papi netcdf) -D_POSIX_C_SOURCE=199309L
CFLAGS=-std=c99 -lm -g -lpapi -fopenmp $(shell pkg-config --cflags libpressio papi netcdf)

ALL=io_test
.PHONY: all clean

io_test: io_test.c
	$(MPICC) $(CFLAGS) -o io_test io_test.c -I$(HDF5_INC) -L$(HDF5_LIB) $(LDLIBS) -lhdf5 -lmpi

szx_io_test: szx_io_test.c
	$(MPICC) $(CFLAGS) -o szx_io_test szx_io_test.c -I$(HDF5_INC) -I$(SZX_INC) -L$(SZX_LIB) -L$(HDF5_LIB) $(LDLIBS) -lhdf5 -lmpi -lSZx

clean:
	rm io_test
