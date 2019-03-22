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

### object files

UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

### programs

vgpvalidate: vgpvalidate.c $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

### end of file
