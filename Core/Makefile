# makefile for vgp-tools/src containing library and core utilities

DEST_DIR = ~/bin

#CFLAGS= -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing -DNDEBUG # NDEBUG drops asserts
CFLAGS= -g -Wall -Wextra -Wno-unused-result -fno-strict-aliasing  # for debugging

LIB = libONE.a
PROGS = ONEstat ONEview

all: $(LIB) $(PROGS)

clean:
	$(RM) *.o ONEstat ONEview $(LIB) ZZ*
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
	$(CC) $(CFLAGS) -o $@ $^

ONEview: ONEview.c utils.o $(LIB)
	$(CC) $(CFLAGS) -o $@ $^

### crude test

TEST:
	./ONEstat -S ../VGP/VGP_1_1.def test.aln
	./ONEview -b -S ../VGP/VGP_1_1.def test.aln > ZZ.1aln
	./ONEstat -S ../VGP/VGP_1_1.def ZZ.1aln
	./ONEstat ZZ.1aln
	./ONEview ZZ.1aln > ZZ.aln
	./ONEstat ZZ.aln
	./ONEstat -C "D C 1 6 STRING" ZZ.aln
	./ONEstat -C "P 3 aln\nD C 1 6 STRING" ZZ.1aln
#	./ONEstat -C "D C 1 3 INT" ZZ.aln
#	./ONEstat -C "D C 2 3 INT 4 CHAR" ZZ.aln
#	./ONEstat -C "D X 1 3 INT" ZZ.aln

### end of file
