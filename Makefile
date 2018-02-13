## Simple makefile to compile BayesTraitsv3 ##
#
# Program info:
# BayesTraits V3.0 (Mar 14 2017)
# Mark Pagel and Andrew Meade
# www.evolution.reading.ac.uk
#
# source: http://www.evolution.rdg.ac.uk/BayesTraitsV3/BayesTraitsV3.html
#
# By default, compiles serial version. To compile threaded version, do:
#
# make threaded

CC = gcc

# *lots* of warnings, so turning off -Wall
#CFLAGS = -O3 -Wall
CFLAGS = -O3
LFLAGS = -lm -lgsl -lgslcblas -lnlopt -llapack

threaded: CFLAGS += -fopenmp -DOPENMP_THR

.PHONY: default all serial threaded clean

default: serial
all: default

OBJS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

serial: $(OBJS)
	$(CC) $(OBJS) -Wall $(LFLAGS) -o BayesTraitsv3_serial

threaded: $(OBJS)
	$(CC) $(OBJS) -Wall $(LFLAGS) -DOPENMP_THR -fopenmp -o BayesTraitsv3_threaded

clean:
	rm -f *.o BayesTraitsv3_*

