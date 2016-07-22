TEST_NEWICK=../../bin/test_newick_tree.py

TREES := $(shell cd trees && find . -name "*.newick" -type f)

.PHONY: all $(TREES)

all: $(TREES)
	echo $(TREES)

$(TREES):
	echo ahoj
	echo $@
	$(TEST_NEWICK) -n ./trees/$@
