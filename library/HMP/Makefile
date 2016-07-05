.PHONY: clean all
.DEFAULT_GOAL := all

include ../conf.mk

all: .complete

.complete.download:
	curl $(HMP_SERVER)/HMREFG/all_seqs.fa.bz2 | bzip2 -d > all_seqs.fa
	touch $@

.complete.split: .complete.download
	mkdir -p ../HMPsplit
	../../bin/split_allseq.py all_seqs.fa ../HMPsplit
	#rm all_seqs.fa
	touch $@

.complete.faidx: .complete.split
	cd ../HMPsplit
	bash -c $(FAIDX)
	cd ../HMP
	touch $@

.complete: .complete.faidx
	touch $@
