CC=gcc
MPICC=mpicc

SZX_DIR := $(HOME)/SZx/build
SZX_INC := $(SZX_DIR)/include
SZX_LIB := $(SZX_DIR)/lib64

HDF5_DIR := /opt/packages/spack/opt/spack/linux-centos8-zen2/gcc-10.2.0/hdf5-1.10.7-s67cr5wg4npwqsrj4ee5c24zh5wklc6l
HDF5_INC := $(HDF5_DIR)/include
HDF5_LIB := $(HDF5_DIR)/lib
LDLIBS=$(shell pkg-config --libs libpressio papi) -D_POSIX_C_SOURCE=199309L
CFLAGS=-std=c99 -lm -g -lpapi $(shell pkg-config --cflags libpressio papi)

ALL=io_test
.PHONY: all clean

io_test: io_test.c
	$(MPICC) $(CFLAGS) -o io_test io_test.c -I$(HDF5_INC) -L$(HDF5_LIB) $(LDLIBS) -lhdf5 

szx_io_test: szx_io_test.c
	$(MPICC) $(CFLAGS) -o szx_io_test szx_io_test.c -I$(HDF5_INC) -I$(SZX_INC) -L$(SZX_LIB) -L$(HDF5_LIB) $(LDLIBS) -lhdf5 -lmpi -lSZx

clean:
	rm io_test
