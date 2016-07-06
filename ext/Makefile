all: time/time

time/time: time/Makefile
	$(MAKE) -C time

time/Makefile:
	(cd time && autoreconf -f -i && ./configure)

