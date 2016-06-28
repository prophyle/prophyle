.PHONY: clean all
.DEFAULT_GOAL := all

all: simulation_bacteria.fq simulation_bacteria.10000.fq simulation_bacteria.100000.fq simulation_bacteria.1000000.fq

simulation_bacteria.fq:
	curl "http://igm.univ-mlv.fr/~brinda/metang/simulation_bacteria.fq.gz" | gunzip -c > simulation_bacteria.fq
	touch .complete

simulation_bacteria.1000.fq: simulation_bacteria.fq
	head -n 4000 simulation_bacteria.fq > simulation_bacteria.10000.fq

simulation_bacteria.10000.fq: simulation_bacteria.fq
	head -n 40000 simulation_bacteria.fq > simulation_bacteria.10000.fq

simulation_bacteria.100000.fq: simulation_bacteria.fq
	head -n 400000 simulation_bacteria.fq > simulation_bacteria.100000.fq

simulation_bacteria.1000000.fq: simulation_bacteria.fq
	head -n 4000000 simulation_bacteria.fq > simulation_bacteria.1000000.fq

clean:
	rm -fr .complete *.fq

