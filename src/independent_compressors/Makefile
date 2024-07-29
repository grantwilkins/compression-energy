CC=gcc
INCLUDE=-I/home/ac.gwilkins/SZx/build/include
LIB=-L/home/ac.gwilkins/SZx/build/lib
CFLAGS=-O3 -std=c99 -Wall -fPIC -shared
LDFLAGS=-lszx -lm
OPENMP=-fopenmp

all: sz_serial.out

szx_serial.out: szx_serial.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LIB) -o $@ $< $(LDFLAGS)

szz_omp.out: szz_omp.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LIB) -o $@ $< $(LDFLAGS) $(OPENMP)

clean:
	rm -f *.out

