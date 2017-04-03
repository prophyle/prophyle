.PHONY: all clean $(TREES)

include ../conf.mk

TREES := $(shell cd $(PROP_DIR)/trees && find . -name "*.nw" -type f)

all: $(TREES)

clean:

$(TREES):
	@echo $@
	$(TEST_NEWICK) -n ./trees/$@
