# makefile for vgp validator, developed on Richard's Mac.

DEST_DIR = ../bin

#CFLAGS= -O3
CFLAGS= -g				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

ALL = vgpvalidate

all: $(ALL)

clean:
	$(RM) *.o vgpvalidate
	$(RM) -r *.dSYM

install:
	cp $(ALL) $(DEST_DIR)

vgprd.o: vgprd.h

### programs

vgpvalidate: vgpvalidate.c vgprd.o
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

### end of file
