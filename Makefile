DEST_DIR = ~/bin

all: README.html
	cd Myers; make; make install
	cd Durbin; make; make install

clean:
	cd Myers; make clean;
	cd Durbin; make clean;
	cd bin; rm *
	rm README.html

%.html: %.md
	pandoc -sS --toc --toc-depth=1 -o $@ $<
