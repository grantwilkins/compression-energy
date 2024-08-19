CC=gcc
SZX_DIR := /jet/home/gwilkins/SZx
INCLUDE=-I/home/ac.gwilkins/SZx/swing_build/include/szx
LIB=-L/home/ac.gwilkins/SZx/swing_build/lib
CFLAGS=-O3 -std=c99 -Wall -fPIC -shared
LDFLAGS=-lm
OPENMP=-fopenmp

all: szx_serial.out szx_omp.out

szx_serial.out: szx_serial.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LIB) -o $@ $< $(LDFLAGS)

szz_omp.out: szz_omp.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LIB) -o $@ $< $(LDFLAGS) $(OPENMP)

clean:
	rm -f *.out

