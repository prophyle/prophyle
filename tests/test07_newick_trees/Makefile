TEST_NEWICK=../../bin/test_newick_tree.py

TREES := $(shell cd trees && find . -name "*.newick" -type f)

.PHONY: all clean $(TREES)

all: $(TREES)

clean:

$(TREES):
	@echo
	@echo $@
	$(TEST_NEWICK) -n ./trees/$@
