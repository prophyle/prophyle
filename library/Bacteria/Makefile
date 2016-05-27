.PHONY: clean all
.DEFAULT_GOAL := all

include ../conf.mk

all: .complete

.complete:
	curl $(FTP_SERVER)/genomes/archive/old_refseq/Bacteria/all.fna.tar.gz | tar xvz
	touch .complete
