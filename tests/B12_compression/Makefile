.PHONY: all clean compr decompr decompr_noklcp

include ../conf.mk

K=2
tree=data/tree.nw
index=_index

all: decompr decompr_noklcp

decompr: compr
	$(PROP) decompress $(index).tar.gz _decompress1
	cmp _decompress1/$(index)/index.fa.bwt $(index)/index.fa.bwt
	cmp _decompress1/$(index)/index.fa.2.klcp $(index)/index.fa.2.klcp

decompr_noklcp: compr
	$(PROP) decompress -K $(index).tar.gz _decompress2
	cmp _decompress2/$(index)/index.fa.bwt $(index)/index.fa.bwt

compr:
	mkdir -p _decompress1 _decompress2
	$(PROP) index -k 2 -A $(tree) $(index)
	$(PROP) compress $(index)

clean:
	rm -rf _*
