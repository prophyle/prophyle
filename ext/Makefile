.PHONY: all clean

all: time/time bwa/bwa kraken-db

time/time: time/Makefile
	$(MAKE) -C time

time/Makefile:
	#(cd time && autoreconf -f -i && ./configure)
	(cd time && ./configure)

bwa/bwa:
	$(MAKE) -C bwa

kraken-db:
	$(MAKE) -C $@

clean:
	make -C time clean
	make -C bwa clean
	make -C kraken-db clean
