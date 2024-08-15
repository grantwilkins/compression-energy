LDLIBS=$(shell pkg-config --libs libpressio)
CFLAGS= -lm -fopenmp -g $(shell pkg-config --cflags libpressio)
ALL=get_compressor_options compress_cpu list_compressors compress_cpu_omp
.PHONY: all clean

all:  $(ALL)

compress_gpu: compress_gpu.c
	nvcc compress_gpu.c $> -g $(shell pkg-config --cflags libpressio) $(LDLIBS) -lm -g

clean:
	$(RM) $(ALL) $(wildcard *.o)
