.PHONY: all clean

include ../conf.mk

all:
	$(PROP) download -d seqs viruses
	$(PROP) index seqs/test_viruses.nw test_index

clean:
	rm -fr test_index seqs
