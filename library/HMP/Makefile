.PHONY: clean all
.DEFAULT_GOAL := all

include ../conf.mk

all: .complete

.complete:
	curl $(HMP_SERVER)/HMREFG/all_seqs.fa.bz2 | bzip2 -d > all_seqs.fa
	curl $(FTP_SERVER)/pub/taxonomy/gi_taxid_nucl.dmp.gz | gzip -d > gi_taxid_nucl.dmp
	touch .complete
