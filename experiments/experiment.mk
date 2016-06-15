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
	TIME?=gtime
else
	TIME?=time
endif
TTIME:=$(TIME) -v

all: index.fa.bwt index.fa.$(K).bit.klcp _time_log.log

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
	$(TTIME) -o 1_kmer_propagation.log \
	$(MAKE) -f Makefile.generated

	touch index/.complete

index.fa: index/.complete
	$(TTIME) -o 2_merging_fasta.log \
	$(FINAL_FA) index > index.fa

%.sa %.pac %.bwt %.amb %.ann: %
	$(TTIME) -o 3_bwa_index.log \
	$(BWA) index -a is $<

%.$(K).bit.klcp: % %.bwt
	$(TTIME) -o 4_klcp.log \
	$(EXK) index -k $(K) $<

%.fai: %
	$(TTIME) -o 5_fasta_index.log \
	$(SAMTOOLS) faidx $<

_time_log.log: index.fa.$(K).bit.klcp
	tail -n +1 [0-9]*.log > _time_log.log

clean:
	rm -f index.fa index.fa.* Makefile.generated *.log
	rm -fr index/

