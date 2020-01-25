DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = VGPseq VGPzip VGPpacbio VGPcloud Dazz2pbr Dazz2sxs VGPpipe

all: deflate.lib $(ALL)

VGP_CORE = ../Durbin/vgprd.c ../Durbin/compression.c
VGP_DPND = $(VGP_CORE) ../Durbin/vgprd.h ../Durbin/compression.h

GENE_CORE = gene_core.c LIBDEFLATE/libdeflate.a
GENE_DPND = $(GENE_CORE) gene_core.h

deflate.lib: LIBDEFLATE
	cd LIBDEFLATE; make; cd ..

VGPzip: VGPzip.c $(GENE_DPND)
	gcc $(CFLAGS) -o VGPzip VGPzip.c $(GENE_CORE) -lpthread

VGPseq: VGPseq.c $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o VGPseq VGPseq.c $(VGP_CORE) $(GENE_CORE) -lpthread

VGPpacbio: VGPpacbio.c pb_expr.c pb_expr.h $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o VGPpacbio VGPpacbio.c pb_expr.c $(VGP_CORE) $(GENE_CORE) -lpthread

VGPpipe: VGPpipe.c $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o VGPpipe VGPpipe.c $(VGP_CORE) $(GENE_CORE) -lpthread

Dazz2pbr: Dazz2pbr.c $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o Dazz2pbr Dazz2pbr.c $(VGP_CORE) $(GENE_CORE) -lpthread

Dazz2sxs: Dazz2sxs.c $(VGP_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o Dazz2sxs Dazz2sxs.c $(VGP_CORE) $(GENE_CORE) -lpthread

VGPcloud: VGPcloud.c radix.c radix.h exsort.c exsort.h $(GENE_DPND)
	gcc $(CFLAGS) -o VGPcloud -DWORD_SIZE=4 VGPcloud.c radix.c exsort.c $(GENE_CORE) -lpthread

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f VGPtools.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf VGPtools.tar.gz Makefile *.h *.c LIBDEFLATE
