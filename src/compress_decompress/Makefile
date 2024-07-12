DATADIR=$(realpath $(shell pwd)/../..)/datasets/
LDLIBS=$(shell pkg-config --libs libpressio)
CFLAGS=-std=c99 -lm -g $(shell pkg-config --cflags libpressio)
CPPFLAGS=-DDATADIR="\"$(DATADIR)\""
ALL=get_compressor_options.out compress_cpu.out list_compressors.out
.PHONY: all clean

all:  $(ALL)
clean:
	$(RM) $(ALL) $(wildcard *.o)