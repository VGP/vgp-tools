DEST_DIR = ~/bin

all: README.html
	cd Core; make
	cd VGP; make

clean:
	cd Core; make clean
	cd VGP; make clean
	rm *.html

%.html: %.md
	pandoc -s --toc --toc-depth=1 -o $@ $<
# 220721 was "pandoc -sS --toc --toc-depth=1 -o $@ $<" but the -S option has gone
