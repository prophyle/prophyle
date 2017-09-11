.PHONY: all test clean prophyle inc

all: prophyle

prophyle:
	$(MAKE) -C prophyle

test:
	$(MAKE) -C prophyle
	$(MAKE) -C tests clean
	$(MAKE) -C tests

clean:
	$(MAKE) -C prophyle clean
	$(MAKE) -C tests clean

inc:
	./prophyle/increment_version.py
