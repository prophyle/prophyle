.PHONY: all clean

all: bwa/bwa kraken-db/.complete

bwa/bwa:
	$(MAKE) -C bwa

kraken-db/.complete:
	$(MAKE) -C kraken-db

clean:
	make -C bwa clean
	make -C kraken-db clean
