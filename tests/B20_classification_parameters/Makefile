.PHONY: all clean classify d L X LX M K h1 h2 c1 c2

include ../conf.mk

K=3
tree=tree.nw
index=index

all: classify

classify: d L X LX M K h1 h2 c1 c2

d: index.complete
	$(PROP) classify $(index) reads.fq | $(SVH) > class.sam

L: index.complete
	$(PROP) classify $(index) reads.fq -L | $(SVH) > class_L.sam

X: index.complete
	$(PROP) classify $(index) reads.fq -X | $(SVH) > class_X.sam

LX: index.complete
	$(PROP) classify $(index) reads.fq -L -X | $(SVH) > class_LX.sam

M: index.complete
	$(PROP) classify $(index) reads.fq -M | $(SVH) > class_M.krak

K: index.complete
	$(PROP) classify $(index) reads.fq -K | $(SVH) > class_K.sam

h1: index.complete
	$(PROP) classify $(index) reads.fq -m h1 | $(SVH) > class_h1.sam

h2: index.complete
	$(PROP) classify $(index) reads.fq -m h2 | $(SVH) > class_h2.sam

c1: index.complete
	$(PROP) classify $(index) reads.fq -m c1 | $(SVH) > class_c1.sam

c2: index.complete
	$(PROP) classify $(index) reads.fq -m c2 | $(SVH) > class_c2.sam

index.complete:
	$(PROP) index -k $(K) -A $(tree) $(index)
	touch $@

clean:
	rm -rf $(index) $(res) $(diff_res) *.complete *.sam *.krak
