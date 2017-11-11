.PHONY: \
	all prophyle clean install hooks \
	test test_repo test_repo_coverage test_parallel test_package \
	pylint coverage \
	inc pypi sha256 \
	docs readme wpypi wconda \
	deppip depconda \
	submodules \
	help

PIP=/usr/bin/env pip
PYTHON=/usr/bin/env python3

ROOT_DIR = $(shell pwd)

###############
# BASIC RULES #
###############

all: prophyle ## Compile ProPhyle

help: ## Print help message
	    @echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

prophyle: ## Compile ProPhyle
prophyle: hooks
	$(MAKE) -C prophyle

clean: ## Clean
	$(PYTHON) setup.py clean --all
	rm -fr _index_test/ _test_*
	$(MAKE) -C prophyle clean
	$(MAKE) -C tests clean
	($(MAKE) -C docs clean || true) > /dev/null 2> /dev/null

install: ## Install ProPhyle using PIP
install: hooks
	$(PIP) uninstall -y prophyle || true
	$(PIP) install .

hooks: ## Install git hooks
	rm -f .git/hooks/*
	@for x in $$(find bin/hooks -type f); do \
		h=$$(basename "$$x"); \
		echo "Installing hook $$h"; \
		y=".git/hooks/$$h"; \
		rm -f "$$y"; \
		ln -s "../../$$x" "$$y"; \
	done

###########
# TESTING #
###########

test: test_repo

coverage: ## Run test coverage analysis
coverage: test_repo_coverage

test_repo: ## Run unit tests & integration from the repo dir
test_repo:
	$(MAKE) -C tests clean
	$(MAKE) -C tests

test_repo_coverage:
	$(MAKE) -C tests clean
	# replace /usr/bin/env/python3 by coverage
	PATH=$$(pwd)/bin/python3_coverage_wrapper:$$PATH $(MAKE) -C tests

test_parallel: ## Run tests in parallel
	$(MAKE) -C tests clean
	$(MAKE) -C tests parallel


test_package: ## Run integration tests from the Python package
	$(MAKE) -C tests clean
	$(MAKE) -C tests B PROP=prophyle

pylint: ## Run PyLint
	$(PYTHON) -m pylint -d prophyle


#############
# RELEASING #
#############

inc: ## Increment version
inc: hooks
	./prophyle/increment_version.py

pypi: ## Upload ProPhyle to PyPI
pypi: hooks
	$(MAKE) clean
	$(PYTHON) setup.py sdist bdist_wheel upload

## Compute sha256 for the PyPI package
sha256:
	s=$$(curl https://pypi.python.org/pypi/prophyle  2>/dev/null| perl -pe 's/#/\n/g' | grep -o 'https.*\.tar\.gz' | xargs curl -L 2>/dev/null | shasum -a 256 | awk '{print $$1;}'); echo $$s; echo $$s | pbcopy


#######################
# DOCUMENTATION & WEB #
#######################

docs: ## Build and open Sphinx documentation
	$(MAKE) -C docs html
	open docs/.build/html/index.html || true

readme: ## Convert README to HTML
readme: hooks
	rst2html.py README.rst > README.html

wconda: ## Open ProPhyle Bioconda webpage
	open https://bioconda.github.io/recipes/prophyle/README.html

wpypi: ## Open ProPhyle PyPI webpage
	open https://pypi.python.org/pypi/prophyle


########################
# INSTALL DEPENDENCIES #
########################

depconda: ## Install dependencies using Conda
	cat requirements.txt | xargs conda install

deppip: ## Install dependencies using PIP
	cat requirements.txt | xargs $(PIP) install


##############
# SUBMODULES #
##############

submodules: ## Download BWA submodule if missing
	$(MAKE) -C prophyle submodules

