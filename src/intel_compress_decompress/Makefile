DATADIR=/ocean/projects/cis240100p/gwilkins/
LDLIBS=$(shell pkg-config --libs libpressio papi)
CFLAGS=-std=c99 -lm -lpapi -fopenmp -g -lpapi $(shell pkg-config --cflags libpressio papi)
CPPFLAGS=-DDATADIR="\"$(DATADIR)\""
ALL=compress_cpu compress_cpu_omp
.PHONY: all clean

all:  $(ALL)

compress_cpu_omp: compress_cpu_omp.c
	mpicc $(CFLAGS) -o $@ $< $(LDLIBS) 

clean:
	$(RM) $(ALL) $(wildcard *.o)
