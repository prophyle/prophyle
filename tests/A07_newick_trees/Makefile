include ../conf.mk
TREES := $(wildcard $(PROP_DIR)/trees/*.nw)

.PHONY: all clean $(TREES)

all: $(TREES)

clean:

$(TREES):
	$(TEST_NEWICK) "$@"
