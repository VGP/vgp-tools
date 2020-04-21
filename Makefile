DEST_DIR = ~/bin

all: README.html
	cd Core; make
	cd VGP; make

clean:
	cd Core; make clean
	cd VGP; make clean
	rm *.html

%.html: %.md
	pandoc -sS --toc --toc-depth=1 -o $@ $<
