.PHONY: all clean classify d L X LX M K h1 h2 c1 c2

include ../conf.mk

K=3
tree=tree.nw
index=index

all: classify

classify: d L X LX M K h1 h2 c1 c2

d: index.complete
	$(PROP) classify $(index) reads.fq > class.sam

L: index.complete
	$(PROP) classify $(index) reads.fq -L > class_L.sam

X: index.complete
	$(PROP) classify $(index) reads.fq -X > class_X.sam

LX: index.complete
	$(PROP) classify $(index) reads.fq -L -X > class_LX.sam

M: index.complete
	$(PROP) classify $(index) reads.fq -M > class_M.krak

K: index.complete
	$(PROP) classify $(index) reads.fq -K > class_K.sam

h1: index.complete
	$(PROP) classify $(index) reads.fq -m h1 > class_h1.sam

h2: index.complete
	$(PROP) classify $(index) reads.fq -m h2 > class_h2.sam

c1: index.complete
	$(PROP) classify $(index) reads.fq -m c1 > class_c1.sam

c2: index.complete
	$(PROP) classify $(index) reads.fq -m c2 > class_c2.sam

index.complete:
	$(PROP) index -k $(K) -A $(tree) $(index)
	touch $@

clean:
	rm -rf $(index) $(res) $(diff_res) *.complete *.sam *.krak
