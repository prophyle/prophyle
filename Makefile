.PHONY: all test clean prophyle inc docs readme

all: prophyle

prophyle:
	$(MAKE) -C prophyle

test:
	$(MAKE) -C prophyle
	$(MAKE) -C tests clean
	$(MAKE) -C tests

pypi:
	/usr/bin/env python3 setup.py sdist bdist_wheel upload

docs:
	$(MAKE) -C docs html
	open docs/.build/html/index.html || true

readme:
	rst2html.py README.rst > README.html

clean:
	rm -fr build dist prophyle.egg-info
	rm -fr _index_test/ _test_*
	$(MAKE) -C prophyle clean
	$(MAKE) -C tests clean
	$(MAKE) -C docs clean

inc:
	./prophyle/increment_version.py
