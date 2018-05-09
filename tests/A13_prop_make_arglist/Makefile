.PHONY: all clean index

include ../conf.mk

K=5

all: index.fa index_noprop.fa
	# with masking
	diff -c index_noprop.fa expected_noprop.txt
	# with propagation
	diff -c index.fa expected.txt

index_noprop.fa:
	$(PROP) index -R -T -k $(K) -g lib tree.nw index_noprop
	$(NORM) -i index_noprop/index.fa > $@

index.fa:
	$(PROP) index -T -k $(K) -g lib tree.nw index
	$(NORM) -i index/index.fa > $@

clean:
	rm -rf index*
