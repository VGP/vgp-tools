DEST_DIR = ../bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = VGPseq VGPzip VGPpacbio VGPcloud Dazz2pbr Dazz2sxs

all: deflate.lib $(ALL)

deflate.lib: LIBDEFLATE
	cd LIBDEFLATE; make; cd ..

VGPzip: VGPzip.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o VGPzip VGPzip.c gene_core.c LIBDEFLATE/libdeflate.a -lm -lpthread -lz

VGPseq: VGPseq.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o VGPseq VGPseq.c gene_core.c LIBDEFLATE/libdeflate.a -lm -lpthread -lz

VGPpacbio: VGPpacbio.c pb_expr.c pb_expr.h gene_core.c gene_core.h
	gcc $(CFLAGS) -o VGPpacbio VGPpacbio.c pb_expr.c gene_core.c LIBDEFLATE/libdeflate.a -lpthread -lz

VGPcloud: VGPcloud.c radix.c radix.h exsort.c exsort.h gene_core.c gene_core.h
	gcc $(CFLAGS) -o VGPcloud -DWORD_SIZE=4 VGPcloud.c radix.c exsort.c gene_core.c -lm -lz -lpthread

Dazz2pbr: Dazz2pbr.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Dazz2pbr Dazz2pbr.c gene_core.c -lm -lz

Dazz2sxs: Dazz2sxs.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Dazz2sxs Dazz2sxs.c gene_core.c -lm -lz

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f VGPtools.tar.gz
	cd LIBDEFLATE; make clean; cd ..

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf VGPtools.tar.gz Makefile *.h *.c LIBDEFLATE
