.PHONY: clean all
.DEFAULT_GOAL := all

include ../conf.mk

all: .complete

.complete.download:
	curl $(FTP_SERVER)/genomes/archive/old_refseq/Plasmids/plasmids.all.fna.tar.gz | tar xvz
	touch $@

.complete.faidx: .complete.download
	bash -c $(FAIDX)
	touch $@

.complete: .complete.faidx
	touch $@
