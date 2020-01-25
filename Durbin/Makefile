# makefile for vgp validator, developed on Richard's Mac.

DEST_DIR = ../bin

CFLAGS= -O3
#CFLAGS= -g				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

ALL = vgpvalidate vgpview

all: $(ALL)

clean:
	$(RM) *.o vgpvalidate vgpview
	$(RM) -r *.dSYM

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf Durbin.Tools.tar.gz *.c *.h Makefile

vgprd.o: vgprd.h

compression.o: compression.h

### programs

vgpvalidate: vgpvalidate.c vgprd.o compression.o
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

vgpview: vgpview.c vgprd.o compression.o
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

### end of file
