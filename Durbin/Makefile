# makefile for modimizer, developed on Richard's Mac.

#CFLAGS= -O3
CFLAGS= -g				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

all: vgpvalidate

clean:
	$(RM) *.o *~ vgpvalidate

vgprd.o: vgprd.h

### programs

vgpvalidate: vgpvalidate.c vgprd.o
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

### end of file
