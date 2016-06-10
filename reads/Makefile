.PHONY: clean all
.DEFAULT_GOAL := all

all: .complete

.complete:
	curl "http://igm.univ-mlv.fr/~brinda/metang/simulation_bacteria.fq.gz" | gunzip -c > simulation_bacteria.fq
	touch .complete

clean:
	rm -fr .complete *.fq

