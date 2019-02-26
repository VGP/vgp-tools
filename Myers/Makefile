DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = VGPseq VGPpacbio VGPcloud Dazz2pbr Dazz2sxs

all: $(ALL) README.html

VGPseq: VGPseq.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o VGPseq VGPseq.c gene_core.c -lm -lz

VGPpacbio: VGPpacbio.c pb_sam.c pb_sam.h pb_expr.c pb_expr.h gene_core.c gene_core.h
	gcc $(CFLAGS) -o VGPpacbio VGPpacbio.c pb_sam.c pb_expr.c gene_core.c -lz

VGPcloud: VGPcloud.c radix.c radix.h gene_core.c gene_core.h
	gcc $(CFLAGS) -o VGPcloud -DWORD_SIZE=4 VGPcloud.c radix.c gene_core.c -lm

Dazz2pbr: Dazz2pbr.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Dazz2pbr Dazz2pbr.c gene_core.c -lm

Dazz2sxs: Dazz2sxs.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Dazz2sxs Dazz2sxs.c gene_core.c -lm -lz

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
