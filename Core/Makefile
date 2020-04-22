# makefile for vgp-tools/src containing library and core utilities

DEST_DIR = ~/bin

#CFLAGS= -O3 -I$(IDIR) -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
CFLAGS= -g -I$(IDIR) -Wall -Wextra -Wno-unused-result -fno-strict-aliasing  # for debugging
#CFLAGS= -03-I$(IDIR) -DOMP -fopenmp	# for OMP parallelisation - doesn't compile on Mac

LIB = libONE.a
PROGS = ONEstat ONEview

all: $(LIB) $(PROGS)

clean:
	$(RM) *.o ONEstat ONEview $(LIB)
	$(RM) -r *.dSYM

install:
	cp $(PROGS) $(DEST_DIR)

package:
	make clean
	tar -zcf ONE-core.tar.gz *.c *.h Makefile

### library

LIB_OBJS = ONElib.o

ONElib.o: ONElib.h

$(LIB): $(LIB_OBJS)
	ar -cr $@ $^
	ranlib $@

### programs

ONEstat: ONEstat.c utils.o $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ -lm -lz

ONEview: ONEview.c utils.o $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ -lm -lz

### end of file
