.PHONY: clean all
.DEFAULT_GOAL := all

include ../conf.mk

all: .complete.1 .complete.2

.complete.1:
	curl $(HMP_SERVER)/HMREFG/all_seqs.fa.bz2 | bzip2 -d > all_seqs.fa
	touch $@

.complete.2: .complete.1
	mkdir ../HMPsplit
	../../bin/split_allseq.py all_seqs.fa ../HMPsplit
	#rm all_seqs.fa
	touch $@
