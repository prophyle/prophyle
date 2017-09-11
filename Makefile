.PHONY: all test clean prophyle inc docs readme

all: prophyle

prophyle:
	$(MAKE) -C prophyle

test:
	$(MAKE) -C prophyle
	$(MAKE) -C tests clean
	$(MAKE) -C tests

docs:
	$(MAKE) -C docs html
	open docs/.build/html/index.html || true

readme:
	rst2html.py README.rst > README.html

clean:
	$(MAKE) -C prophyle clean
	$(MAKE) -C tests clean
	$(MAKE) -C docs clean

inc:
	./prophyle/increment_version.py
