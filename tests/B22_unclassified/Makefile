.PHONY: all clean

include ../conf.mk

K=3
tree=tree.nw
index=index

all: output_wo_header.sam
	test -s $<

output_wo_header.sam: index.complete
	$(PROP) classify $(index) to_be_unclassified.fa | samtools view > $@

index.complete:
	$(PROP) index -k $(K) -A $(tree) $(index)
	touch $@

clean:
	rm -rf $(index)
	rm -f *.complete *.sam

