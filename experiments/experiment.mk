.PHONY: all renderfasta

include *.mk

#ASSEMBLER=../../bin/assembler
EXK=../../bin/exk
BUILD_FA=../../src/build_index.py
BWA=../../bin/bwa
SAMTOOLS=samtools
NEWICK2MAKEFILE=../../bin/newick2makefile.py
FINAL_FA=../../bin/create_final_fasta.py

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	TIME=gtime
else
	TIME=time
endif

all: index.fa.bwt index.fa.$(K).bit.klcp

index/:
	mkdir index

Makefile.generated: $(NEWICK2MAKEFILE) $(TREE)
	$(NEWICK2MAKEFILE) \
	-n $(TREE) \
	-o ./index \
	-l ../../ \
	-k $(K) \
	> Makefile.generated


index/.complete: index/ Makefile.generated
	$(TIME) -o 1_kmer_propagation.log \
	$(MAKE) -f Makefile.generated

	touch index/.complete

index.fa: index/.complete
	$(TIME) -o 2_merging_fasta.log \
	$(FINAL_FA) index > index.fa

	#@cat index/*.reduced.fa > index.fa

%.sa %.pac %.bwt %.amb %.ann: %
	$(TIME) -o 3_bwa_index.log \
	$(BWA) index -a is $<

%.$(K).bit.klcp: % %.bwt
	$(TIME) -o 4_klcp.log \
	$(EXK) index -k $(K) $<

%.fai: %
	$(TIME) -o 5_fasta_index.log \
	$(SAMTOOLS) faidx $<

clean:
	rm -f index.fa index.fa.* Makefile.generated *.log
	rm -fr index/
