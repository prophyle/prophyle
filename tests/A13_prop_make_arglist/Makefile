.PHONY: all clean index

include ../conf.mk

K=5

all: index index_mask
	# with masking
	diff -c index_mask/index.fa expected_masked.txt
	# w/o masking
	diff -c index/index.fa expected.txt

index_mask:
	$(PROP) index -M -T -k $(K) tree.nw $@

index:
	$(PROP) index -T -k $(K) tree.nw $@

clean:
	rm -rf index index_mask
