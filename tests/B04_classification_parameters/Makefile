.PHONY: all clean classify

include ../conf.mk

K=3
tree=tree.nw
index=index

all: classify

classify: index.complete
	$(PROP) classify $(index) reads.fq > class.sam
	$(PROP) classify $(index) reads.fq -L > class_L.sam
	$(PROP) classify $(index) reads.fq -X > class_X.sam
	$(PROP) classify $(index) reads.fq -L -X > class_LX.sam
	$(PROP) classify $(index) reads.fq -M > class_M.krak
	$(PROP) classify $(index) reads.fq -K > class_K.sam
	$(PROP) classify $(index) reads.fq -m h1 > class_h1.sam
	$(PROP) classify $(index) reads.fq -m h2 > class_h2.sam
	$(PROP) classify $(index) reads.fq -m c1 > class_c1.sam
	$(PROP) classify $(index) reads.fq -m c2 > class_c2.sam

index.complete:
	$(PROP) index -k $(K) -A $(tree) $(index)
	touch $@

clean:
	rm -rf $(index) $(res) $(diff_res) *.complete *.sam
