CC=gcc
CXX=g++
MPICC=mpicc
MPICXX=mpicxx

SZX_DIR := /jet/home/gwilkins/SZx
SZX_INC := $(SZX_DIR)/include
SZX_LIB := $(SZX_DIR)/lib

LIBPRESSIO_CFLAGS := $(shell pkg-config --cflagslibpressio)
LIBPRESSIO_LDLIBS := $(shell pkg-config --libs libpressio)

MGARDX_DIR := /jet/home/gwilkins/MGARD/build-openmp-cpu/mgard
MGARDX_INC := $(MGARDX_DIR)/include
MGARDX_LIB := $(MGARDX_DIR)/lib64

DATADIR=/ocean/projects/cis240100p/gwilkins/
LDLIBS=$(shell pkg-config --libs papi)
CFLAGS=-std=c99 -lm -g -lpapi $(shell pkg-config --cflags papi)
CPPFLAGS=-DDATADIR="\"$(DATADIR)\""
ALL=compress_cpu compress_cpu_omp szx_omp szx_serial mgard-x mgard-x_omp

.PHONY: all clean

all: $(ALL)

compress_cpu: compress_cpu.c
	$(CC) $(CFLAGS) $(LIBPRESSIO_CFLAGS) -o compress_cpu compress_cpu.c $(LDLIBS) $(LIBPRESSIO_LDLIBS)

compress_cpu_omp: compress_cpu_omp.c
	$(MPICC) $(CFLAGS) $(LIBPRESSIO_CFLAGS)  -o compress_cpu_omp compress_cpu_omp.c $(LDLIBS) $(LIBPRESSIO_LDLIBS) -fopenmp

szx_omp: szx_omp.c
	$(MPICC) $(CFLAGS) -o szx_omp szx_omp.c -I$(SZX_INC) -L$(SZX_LIB) $(LDLIBS) -lSZx -fopenmp $(LDLIBS) 

szx_serial: szx_serial.c
	$(CC) $(CFLAGS) -o szx_serial szx_serial.c -I$(SZX_INC) -L$(SZX_LIB) $(LDLIBS) -lSZx

mgardx: mgardx.cpp
	$(CXX) $(CFLAGS) -o mgardx mgardx.cpp -I$(MGARDX_INC) -L$(MGARDX_LIB) $(LDLIBS) -lmgard

mgardx_omp: mgardx-omp.cpp
	$(MPICXX) $(CFLAGS) -o mgardx-omp mgardx-omp.cpp -I$(MGARDX_INC) -L$(MGARDX_LIB) $(LDLIBS) -fopenmp -lmgard

clean:
	rm -f compress_cpu compress_cpu_omp szx_omp szx_serial mgardx mgard-x-omp
	rm -f *.o
