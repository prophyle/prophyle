.PHONY: all clean sim abund diff

include ../conf.mk

K=2
tree=data/tree.nw
index=_index
exp=expected_res

all: sim abund diff

abund: sim
	$(PROP) profile _index data/n1.bam _abund -a 1 -l 1
	diff -c _abund $(exp)/abundances_n1.tsv

sim:
	$(PROP) index -k $(K) -A $(tree) $(index)
	$(PROP) similarity $(index) -g data -l 2 -x 0 -n 10 -s artillumina -P -R
	cmp $(index)/sim_mat.npy $(exp)/sim_mat.npy

clean:
	rm -rf _*
	rm -rf data/*.fq data/*.bam
	find -H data -maxdepth 1 -mindepth 1 -type d -print0 | xargs -n1 -0 rm -rf
