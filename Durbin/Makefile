# makefile for Richard's VGP code

DEST_DIR = ~/bin

#CFLAGS= -O3
CFLAGS= -g				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

ALL = 

all: $(ALL)

clean:
	$(RM) *.o $(ALL)
	$(RM) -r *.dSYM

install:
#	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf Durbin.Tools.tar.gz *.c *.h Makefile

### programs

### end of file
