.PHONY: all clean index

include ../conf.mk

K=5

all: index.fa index_masked.fa
	# with masking
	diff -c index_masked.fa expected_masked.txt
	# w/o masking
	diff -c index.fa expected.txt

index_masked.fa:
	$(PROP) index -M -T -k $(K) -g lib tree.nw index_masked && $(NORM) -i index_masked/index.fa > $@

index.fa:
	$(PROP) index -T -k $(K) -g lib tree.nw index && $(NORM) -i index/index.fa > $@

clean:
	rm -rf index*
