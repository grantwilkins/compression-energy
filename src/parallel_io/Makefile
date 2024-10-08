CC=gcc
MPICC=mpicc

SZX_DIR := /home1/10191/gfw/SZx/build
SZX_INC := $(SZX_DIR)/include
SZX_LIB := $(SZX_DIR)/lib

HDF5_INC=$(TACC_HDF5_INC)
HDF5_LIB=$(TACC_HDF5_LIB)
LDLIBS=$(shell pkg-config --libs libpressio papi netcdf) -D_POSIX_C_SOURCE=199309L
CFLAGS=-std=c99 -lm -g -lpapi -fopenmp $(shell pkg-config --cflags libpressio papi netcdf)

ALL=parallel_io
.PHONY: all clean

parallel_io: parallel_io.c
	$(MPICC) $(CFLAGS) -o parallel_io parallel_io.c -I$(HDF5_INC) -L$(HDF5_LIB) $(LDLIBS) -lhdf5 -lmpi

clean:
	rm parallel_io
