.PHONY: all clean

all: time/time bwa/bwa kraken-db/.complete

time/time: time/Makefile
	$(MAKE) -C time

time/Makefile:
	#(cd time && autoreconf -f -i && ./configure)
	(cd time && ./configure)

bwa/bwa:
	$(MAKE) -C bwa

kraken-db/.complete:
	$(MAKE) -C kraken-db

clean:
	make -C time clean
	make -C bwa clean
	make -C kraken-db clean
