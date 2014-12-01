INCLUDES = -I../htslib
CC = gcc
OPTS = -Wall -O2 -g

.PHONY:	all htslib clean clean-all

.SUFFIXES:	.c .o 

.c.o:
	$(CC) -c $(OPTS) $(INCLUDES) $< -o $@

all: bam2wig

htslib:
	$(MAKE) -C ../htslib

bam2wig: htslib bam2wig.o
	$(CC) $(OPTS) $(INCLUDES) -o bam2wig bam2wig.o ../htslib/libhts.a -lz -lpthread

clean:
	rm -f *.o bam2wig

clean-all: clean
	make --directory=../htslib clean
