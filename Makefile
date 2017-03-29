.PHONY: all test

VERBOSE=1

include prophyle/get_nb_jobs.mk

all:
	make -C prophyle
	make -C library
	make -C reads

test:
	make -C prophyle
	make -C tests clean
	make -C tests
