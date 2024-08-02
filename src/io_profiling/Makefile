LDLIBS=$(shell pkg-config --libs libpressio)
CFLAGS=-std=c99 -lm -fopenmp -g $(shell pkg-config --cflags libpressio)
ALL=io_test
.PHONY: all clean

all:  $(ALL)
clean:
	$(RM) $(ALL) $(wildcard *.o)
