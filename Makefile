.PHONY: \
	all prophyle clean install \
	test test_repo test_package \
	inc pypi \
	docs readme wpypi wconda


###############
# BASIC RULES #
###############

all: prophyle

prophyle:
	$(MAKE) -C prophyle

clean:
	rm -fr build dist prophyle.egg-info
	rm -fr _index_test/ _test_*
	$(MAKE) -C prophyle clean
	$(MAKE) -C tests clean
	($(MAKE) -C docs clean || True) > /dev/null 2> /dev/null

install:
	pip install --upgrade .


###########
# TESTING #
###########

test: test_repo

# unit tests & integration, invoked in the repo dir
test_repo:
	$(MAKE) -C prophyle
	$(MAKE) -C tests clean
	$(MAKE) -C tests

# integration tests, invoked from the pip package dir
test_package:
	$(MAKE) -C tests clean
	$(MAKE) -C tests B PROP=prophyle


#############
# RELEASING #
#############

inc:
	./prophyle/increment_version.py

pypi:
	make clean
	/usr/bin/env python3 setup.py sdist bdist_wheel upload


#######################
# DOCUMENTATION & WEB #
#######################

docs:
	$(MAKE) -C docs html
	open docs/.build/html/index.html || true

readme:
	rst2html.py README.rst > README.html

wconda:
	open https://bioconda.github.io/recipes/prophyle/README.html

wpypi:
	open https://pypi.python.org/pypi/prophyle

