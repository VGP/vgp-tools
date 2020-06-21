DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = VGPzip VGPseq VGPpair VGPpacbio VGPcloud Dazz2pbr Dazz2sxs

all: deflate.lib libhts.a $(ALL)

ONE_LIB = ../Core/ONElib.c
ONE_DPND = $(ONE_LIB) ../Core/ONElib.h VGPschema.h

GENE_CORE = gene_core.c LIBDEFLATE/libdeflate.a
GENE_DPND = $(GENE_CORE) gene_core.h

deflate.lib: LIBDEFLATE
	cd LIBDEFLATE; make; cd ..

libhts.a: HTSLIB
	cd HTSLIB; make; cd ..

VGPzip: VGPzip.c $(GENE_DPND)
	gcc $(CFLAGS) -o VGPzip VGPzip.c $(GENE_CORE) -lpthread

VGPseq: VGPseq.c $(ONE_DPND) $(GENE_DPND) HTSLIB/libhts.a
	gcc $(CFLAGS) -o VGPseq -I./HTSLIB VGPseq.c $(ONE_LIB) $(GENE_CORE) HTSLIB/libhts.a -lpthread -lz -lm -lbz2 -llzma -lcurl

VGPpair: VGPpair.c $(ONE_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o VGPpair VGPpair.c $(ONE_LIB) $(GENE_CORE) -lpthread

VGPpacbio: VGPpacbio.c pb_expr.c pb_expr.h $(ONE_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o VGPpacbio VGPpacbio.c pb_expr.c $(ONE_LIB) $(GENE_CORE) -lpthread

VGPcloud: VGPcloud.c msd.sort.c msd.sort.h $(ONE_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o VGPcloud VGPcloud.c msd.sort.c $(ONE_LIB) $(GENE_CORE) -lpthread

Dazz2pbr: Dazz2pbr.c $(ONE_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o Dazz2pbr Dazz2pbr.c $(ONE_LIB) $(GENE_CORE) -lpthread

Dazz2sxs: Dazz2sxs.c $(ONE_DPND) $(GENE_DPND)
	gcc $(CFLAGS) -o Dazz2sxs Dazz2sxs.c $(ONE_LIB) $(GENE_CORE) -lpthread

tidyup:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f VGPtools.tar.gz

clean:
	cd LIBDEFLATE; make clean; cd ..
	cd HTSLIB; make clean; cd ..
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f VGPtools.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf VGPtools.tar.gz Makefile *.h *.c LIBDEFLATE HTSLIB
