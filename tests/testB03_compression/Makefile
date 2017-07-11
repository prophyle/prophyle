.PHONY: all clean

include ../conf.mk

K=2
tree=tree.nw
index=index

all:
	$(PROP) index -k 2 $(tree) $(index)
	$(PROP) compress $(index)
	mkdir -p decompress
	$(PROP) decompress $(index).tar.gz decompress
	cmp decompress/$(index)/index.fa.bwt $(index)/index.fa.bwt


clean:
	rm -rf $(index) $(res) $(diff_res)
