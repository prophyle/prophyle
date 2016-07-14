.PHONY: all test

VERBOSE=1

include src/get_nb_jobs.mk

all:
	make -C src
	make -C library
	make -C reads
	./src/run_serial.sh -j $(JOBS)

test:
	make -C src NOTIME=1
	make -C tests clean
	make -C tests