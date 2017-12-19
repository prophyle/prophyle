# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = ProPhyle
SOURCEDIR     = .
BUILDDIR      = .build

WEBSITEDIR    = ../../prophyle.github.io

all: html

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile clean sync view

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).

sync:
	$(MAKE) html
	rm -fr $(WEBSITEDIR)/*
	cp -va .build/html/* $(WEBSITEDIR)
	cd $(WEBSITEDIR) && git add '*' && git diff HEAD && (git commit || true) && git push 

clean:
	rm -fr .build

view:
	open .build/html/index.html

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
