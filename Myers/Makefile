DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = VGPzip VGPseq VGPpair VGPpacbio VGPcloud Dazz2pbr Dazz2sxs VGPpipe

all: deflate.lib $(ALL)

VGP_LIB = -L../lib -lVGP
VGP_DPND = ../include/VGPlib.h

GENE_CORE = gene_core.c LIBDEFLATE/libdeflate.a
GENE_DPND = $(GENE_CORE) gene_core.h

deflate.lib: LIBDEFLATE
	cd LIBDEFLATE; make; cd ..

VGPzip: VGPzip.c $(GENE_DPND)
	gcc $(CFLAGS) -o VGPzip VGPzip.c $(GENE_CORE) -lpthread

VGPseq: VGPseq.c $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o VGPseq VGPseq.c $(VGP_LIB) $(GENE_CORE) -lpthread

VGPpair: VGPpair.c $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o VGPpair VGPpair.c $(VGP_LIB) $(GENE_CORE) -lpthread

VGPpacbio: VGPpacbio.c pb_expr.c pb_expr.h $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o VGPpacbio VGPpacbio.c pb_expr.c $(VGP_LIB) $(GENE_CORE) -lpthread

VGPpipe: VGPpipe.c $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o VGPpipe VGPpipe.c $(VGP_LIB) $(GENE_CORE) -lpthread

Dazz2pbr: Dazz2pbr.c $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o Dazz2pbr Dazz2pbr.c $(VGP_LIB) $(GENE_CORE) -lpthread

Dazz2sxs: Dazz2sxs.c $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o Dazz2sxs Dazz2sxs.c $(VGP_LIB) $(GENE_CORE) -lpthread

VGPcloud: VGPcloud.c lsd.sort.c lsd.sort.h exsort.c exsort.h $(GENE_DPND)
	gcc $(CFLAGS) -o VGPcloud VGPcloud.c lsd.sort.c exsort.c $(GENE_CORE) -lpthread

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f VGPtools.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf VGPtools.tar.gz Makefile *.h *.c LIBDEFLATE
