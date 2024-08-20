CC=gcc
CXX=g++
MPICC=mpicc
MPICXX=mpicxx

SZX_DIR := /jet/home/gwilkins/SZx
SZX_INC := $(SZX_DIR)/include
SZX_LIB := $(SZX_DIR)/lib

DATADIR=/ocean/projects/cis240100p/gwilkins/
LDLIBS=$(shell pkg-config --libs libpressio papi)
CFLAGS=-std=c99 -lm -g -lpapi -D_POSIX_C_SOURCE=199309L $(shell pkg-config --cflags libpressio papi)
CPPFLAGS=-DDATADIR="\"$(DATADIR)\""
ALL=compress_cpu compress_cpu_omp szx_omp szx_serial

.PHONY: all clean

all: $(ALL)

compress_cpu: compress_cpu.c
	$(CC) $(CFLAGS) -o compress_cpu compress_cpu.c $(LDLIBS)

compress_cpu_omp: compress_cpu_omp.c
	$(MPICC) $(CFLAGS) -o compress_cpu_omp compress_cpu_omp.c $(LDLIBS) -fopenmp

szx_omp: szx_omp.c
	$(MPICC) $(CFLAGS) -o szx_omp szx_omp.c -I$(SZX_INC) -L$(SZX_LIB) $(LDLIBS) -lSZx -fopenmp

szx_serial: szx_serial.c
	$(CC) $(CFLAGS) -o szx_serial szx_serial.c -I$(SZX_INC) -L$(SZX_LIB) -lSZx $(LDLIBS) 
clean:
	rm -f compress_cpu compress_cpu_omp
	rm -f *.o
