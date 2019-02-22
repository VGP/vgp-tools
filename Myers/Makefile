PATH_HDF5 = /sw
PATH_HDF5 = /usr/local/hdf5

DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = VGPpair VGPpacbio VGPcloud

all: $(ALL) README.html

VGPpair: VGPpair.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o VGPpair VGPpair.c compression.c gene_core.c -lm -lz

VGPpacbio: VGPpacbio.c pb_sam.c pb_sam.h pb_expr.c pb_expr.h gene_core.c gene_core.h
	gcc $(CFLAGS) -I$(PATH_HDF5)/include -L$(PATH_HDF5)/lib -o VGPpacbio VGPpacbio.c pb_sam.c pb_expr.c gene_core.c -lhdf5 -lz

VGPcloud: VGPcloud.c radix.c radix.h gene_core.c gene_core.h
	gcc $(CFLAGS) -o VGPcloud -DWORD_SIZE=4 VGPcloud.c radix.c gene_core.c -lm

README.html: ../README.md
	pandoc -sS --toc --toc-depth=1 -o README.html <../README.md

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f VGPtools.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf VGPtools.tar.gz README.md Makefile *.h *.c
