# makefile for vgp-tools/src containing library and core utilities

DEST_DIR = ~/bin

IDIR = ../include

#CFLAGS= -O3 -I$(IDIR) -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
CFLAGS= -g -I$(IDIR) -Wall -Wextra -Wno-unused-result -fno-strict-aliasing  # for debugging
#CFLAGS= -03-I$(IDIR) -DOMP -fopenmp	# for OMP parallelisation - doesn't compile on Mac

LIB = ../lib/libVGP.a
PROGS = VGPstat VGPview

all: $(LIB) $(PROGS)

clean:
	$(RM) *.o VGPstat VGPview
	$(RM) -r *.dSYM

install:
	cp $(PROGS) $(DEST_DIR)

package:
	make clean
	tar -zcf VGP-core.tar.gz *.c *.h Makefile

### library

LIB_OBJS = VGPlib.o compression.o

VGPlib.o: $(IDIR)/VGPlib.h
compression.o: compression.h

$(LIB): $(LIB_OBJS)
	ar -cr $@ $^
	ranlib $@

### programs

VGPstat: VGPstat.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ -lm -lz

VGPview: VGPview.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ -lm -lz

### end of file
