LDLIBS=$(shell pkg-config --libs libpressio) -lzfp
CFLAGS= -lm -fopenmp -g $(shell pkg-config --cflags libpressio)
ALL=get_compressor_options compress_cpu list_compressors compress_cpu_omp compress_gpu

.PHONY: all clean

all: $(ALL)

compress_gpu: compress_gpu.c
    nvcc compress_gpu.c -o $@ -g $(shell pkg-config --cflags libpressio) $(LDLIBS) -lm -g -lnvml

clean:
    $(RM) $(ALL) $(wildcard *.o)