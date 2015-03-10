INCLUDES = -I../htslib
CC = gcc
OPTS = -Wall -O3

.PHONY: all htslib clean clean-all

.SUFFIXES: .c .o

.c.o:
	$(CC) -c $(OPTS) $(INCLUDES) $< -o $@

all: ConvertPhredQuals

htslib:
	$(MAKE) -C ../htslib

ConvertPhredQuals: htslib ConvertPhredQuals.o
	$(CC) $(OPTS) -o ConvertPhredQuals ConvertPhredQuals.o ../htslib/libhts.a -lz -lpthread

clean:
	rm -f *.o ConvertPhredQuals

clean-all: clean
	make --directory=../htslib clean
