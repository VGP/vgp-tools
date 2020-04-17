DEST_DIR = ~/bin

all: README.html
	cd src; make; make install
	cd Myers; make; make install
	cd Durbin; make; make install

clean:
	cd src; make clean
	cd Myers; make clean
	cd Durbin; make clean
	rm *.html

%.html: %.md
	pandoc -sS --toc --toc-depth=1 -o $@ $<
