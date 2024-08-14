DATADIR=/ocean/projects/cis240100p/gwilkins/
LDLIBS=$(shell pkg-config --libs libpressio)
CFLAGS=-std=c99 -lm -lpapi -fopenmp -g -lpapi $(shell pkg-config --cflags libpressio)
CPPFLAGS=-DDATADIR="\"$(DATADIR)\""
ALL=get_compressor_options compress_cpu list_compressors compress_cpu_omp
.PHONY: all clean

all:  $(ALL)
clean:
	$(RM) $(ALL) $(wildcard *.o)
