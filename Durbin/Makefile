# makefile for vgp validator, developed on Richard's Mac.

DEST_DIR = ../bin

#CFLAGS= -O3
CFLAGS= -g				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

ALL = vgpvalidate vgpview

all: $(ALL)

clean:
	$(RM) *.o vgpvalidate
	$(RM) -r *.dSYM

install:
	cp $(ALL) $(DEST_DIR)

vgprd.o: vgprd.h

compression.o: ../Myers/compression.h
	$(CC) $(CFLAGS) -c ../Myers/compression.c -o ./compression.o

### programs

vgpvalidate: vgpvalidate.c vgprd.o compression.o
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

vgpview: vgpview.c vgprd.o compression.o
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

### end of file
