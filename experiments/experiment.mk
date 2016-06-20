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
#TTIME:=$(TIME) -v
TTIME:=$(TIME)

all: index.fa.bwt index.fa.$(K).bit.klcp _main_log.log

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
	$(BWA) index $<

%.$(K).bit.klcp: % %.bwt
	$(TTIME) -o 4_klcp.log \
	$(EXK) index -k $(K) $<

#%.fai: %
#	$(TTIME) -o x_fasta_index.log \
#	$(SAMTOOLS) faidx $<

kmers_rolling.txt: ../../reads/simulation_bacteria.1000000.fq index.fa.bwt
	$(TTIME) -o 5_matching_rolling.log \
	$(EXK) match -k $(K) -u -v index.fa ../../reads/simulation_bacteria.1000000.fq > kmers_rolling.txt

kmers_restarted.txt: ../../reads/simulation_bacteria.1000000.fq index.fa.bwt kmers_rolling.txt
	$(TTIME) -o 6_matching_restarted.log \
	$(EXK) match -k $(K) -v index.fa ../../reads/simulation_bacteria.1000000.fq > kmers_restarted.txt

_main_log.log: index.fa.$(K).bit.klcp kmers_rolling.txt kmers_restarted.txt
	du -sh *.fa.* | grep -v "fa.amb" > 7_index_size.log
	echo > _main_log.log
	date >> _main_log.log
	pwd >> _main_log.log
	echo >> _main_log.log

	tail -n +1 [0-9]*.log >> _main_log.log

clean:
	rm -f index.fa index.fa.* Makefile.generated *.log
	rm -fr index/

