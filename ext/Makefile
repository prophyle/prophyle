.PHONY: all clean

all: time/time bwa/bwa

time/time: time/Makefile
	$(MAKE) -C time

time/Makefile:
	#(cd time && autoreconf -f -i && ./configure)
	(cd time && ./configure)

bwa/bwa:
	$(MAKE) -C bwa

clean:
	make -C time clean
	make -C bwa clean
