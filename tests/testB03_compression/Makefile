.PHONY: all clean

include ../conf.mk

K=2
tree=tree.nw
index=index

all:
	$(PROP) index -k 2 $(tree) $(index)
	$(PROP) compress $(index)

	mkdir -p decompress1 decompress2

	$(PROP) decompress $(index).tar.gz decompress1
	cmp decompress1/$(index)/index.fa.bwt $(index)/index.fa.bwt

	$(PROP) decompress -K $(index).tar.gz decompress2
	cmp decompress2/$(index)/index.fa.bwt $(index)/index.fa.bwt

clean:
	rm -rf $(index) $(res) $(diff_res)
