CC=gcc
CXX=g++
MPICC=mpicc
MPICXX=mpicxx

SZX_DIR := /home1/10191/gfw/SZx/build
SZX_INC := $(SZX_DIR)/include
SZX_LIB := $(SZX_DIR)/lib

CUSZP_DIR := /jet/home/gwilkins/cuSZp
CUSZP_INC := $(CUSZP_DIR)/include
CUSZP_LIB := $(CUSZP_DIR)/build


LIBPRESSIO_CFLAGS := $(shell pkg-config --cflags libpressio)
LIBPRESSIO_LDLIBS := $(shell pkg-config --libs libpressio)

LDLIBS=$(shell pkg-config --libs papi) -D_POSIX_C_SOURCE=199309L
CFLAGS= -lm -g -lpapi $(shell pkg-config --cflags papi)
ALL=compress_cpu compress_cpu_omp szx_omp szx_serial
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

cuszp: cuszp.cu
	nvcc $(CFLAGS) -o cuszp cuszp.cu -I$(CUSZP_INC) -L$(CUSZP_LIB) $(LDLIBS) -lcuSZp -lnvidia-ml

clean:
	rm -f compress_cpu compress_cpu_omp szx_omp szx_serial mgardx mgard-x-omp
	rm -f *.o
