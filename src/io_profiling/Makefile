LDLIBS=$(shell pkg-config --libs libpressio papi)
CFLAGS=-std=c99 -lm -fopenmp -g -lpapi $(shell pkg-config --cflags libpressio papi)
ALL=io_test
.PHONY: all clean

all:  $(ALL)
clean:
	$(RM) $(ALL) $(wildcard *.o)
