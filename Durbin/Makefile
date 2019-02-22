# makefile for modimizer, developed on Richard's Mac.

#CFLAGS= -O3
CFLAGS= -g				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

all: vgpvalidate

clean:
	$(RM) *.o *~ vgpvalidate

### object files

UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

### programs

vgpvalidate: vgpvalidate.c $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

### end of file
