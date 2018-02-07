.PHONY: all clean autocomplete_full autocomplete_partial

include ../conf.mk

K=2
tree1=data/tree.nw
tree2=tree_simple.nw
index1=_index1
index2=_index2

all: autocomplete_full autocomplete_partial

autocomplete_partial:
	$(PROP) index -A -k 2 $(tree1) $(index1)

autocomplete_full:
	$(PROP) index -A -k 2 -g data $(tree2) $(index2)

clean:
	rm -rf _*
