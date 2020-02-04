# makefile for vgp validator, developed on Richard's Mac.

DEST_DIR = ../bin

#CFLAGS= -O3
CFLAGS= -g				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

ALL = VGPstat VGPview

all: $(ALL)

clean:
	$(RM) *.o VGPstat VGPview
	$(RM) -r *.dSYM

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf Durbin.Tools.tar.gz *.c *.h Makefile

VGPlib.o: compression.c VGPlib.h compression.h

### programs

VGPstat: VGPstat.c VGPlib.o
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

VGPview: VGPview.c VGPlib.o
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

### end of file
