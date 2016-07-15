.PHONY: all renderfasta

V=1
SHELL:=/bin/bash

include *.mk

#ASSEMBLER=../../bin/assembler
EXK=../../bin/exk
BUILD_FA=../../src/build_index.py
BWA=../../bin/bwa
SAMTOOLS?=samtools
NEWICK2MAKEFILE=../../bin/newick2makefile.py
ASSIGNMENT=../../bin/read_assignment.py
FINAL_FA=../../bin/create_final_fasta.py

READS?=../../reads/simulation_bacteria.1000000.fq
KLCP=index.fa.$(K).bit.klcp
#UNAME_S := $(shell uname -s)
#
#ifeq ($(UNAME_S),Darwin)
#	TIME?=gtime
#else
#	TIME?=time
#endif
#TTIME:=$(TIME) -v
#TTIME:=$(TIME)

include ../../src/get_nb_jobs.mk

TIME=../../bin/time
TTIME:=DATETIME=`date` && $(TIME) -f "$${DATETIME}\njobs: $(JOBS)\n%C\n%Uuser %Ssystem %Eelapsed %PCPU (%Xavgtext+%Davgdata %Mmaxresident)k\n%Iinputs+%Ooutputs (%Fmajor+%Rminor)pagefaults %Wswaps"

ifdef MASK_REPEATS
	REP_PARAM=-r
else
	REP_PARAM=
endif



all: index.fa.sa index.fa.$(K).bit.klcp _main_log.log _main_log.md \
	assigned_reads.txt assigned_reads_simlca.txt

index/.complete: $(TREE)
	mkdir -p index

	$(NEWICK2MAKEFILE) \
	-n $(TREE) \
	-o ./index \
	-l ../../ \
	-k $(K) \
	$(REP_PARAM) \
	> Makefile.generated

	$(TTIME) -o 1.1_kmer_propagation.log \
	$(MAKE) -f Makefile.generated

	touch index/.complete

index.fa: index/.complete
	$(TTIME) -o 1.2_merging_fasta.log \
	$(FINAL_FA) index > index.fa

index.fa.pac: index.fa 
	$(TTIME) -o 2.1_bwa_fa2pac.log \
	$(BWA) fa2pac index.fa index.fa

index.fa.bwt: index.fa.pac 
	$(TTIME) -o 2.2_bwa_pac2bwtgen.log \
	$(BWA) pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

	$(TTIME) -o 2.3_bwa_bwtupdate.log \
	$(BWA) bwtupdate index.fa.bwt

index.fa.sa: index.fa.bwt 
	$(TTIME) -o 2.4_bwa_bwt2sa.log \
	$(BWA) bwt2sa index.fa.bwt index.fa.sa

$(KLCP): index.fa index.fa.bwt index.fa.sa
	$(TTIME) -o 2.5_klcp.log \
	$(EXK) index -k $(K) $<

%.fai: %
	$(SAMTOOLS) faidx $<

kmers_rolling.txt: index.fa.sa $(KLCP) 
	$(TTIME) -o 3.1a_matching_rolling.log \
	$(EXK) match -l 3.1b_matching_rolling.log  \
		-k $(K) -u index.fa $(READS) > $@

kmers_restarted.txt: $(READS) index.fa.sa \
	kmers_rolling.txt
	$(TTIME) -o 3.2a_matching_restarted.log \
	$(EXK) match -l 3.2b_matching_restarted.log \
		-k $(K) index.fa $(READS) > $@

kmers_rolling_skipping.txt: $(READS) index.fa.sa $(KLCP) \
	kmers_restarted.txt
	$(TTIME) -o 3.3a_matching_rolling_skipping.log \
	$(EXK) match -l 3.3b_matching_rolling_skipping.log \
		-k $(K) -u -s index.fa $(READS) > $@

kmers_restarted_skipping.txt: $(READS) index.fa.sa kmers_rolling.txt \
	kmers_rolling_skipping.txt
	$(TTIME) -o 3.4a_matching_restarted_skipping.log \
	$(EXK) match -l 3.4b_matching_restarted_skipping.log \
		-k $(K) -s index.fa $(READS) > $@

assigned_reads.txt: kmers_rolling.txt $(TREE)
	$(TTIME) -o 4.1_read_assignment.log \
	$(ASSIGNMENT) -i $< -n $(TREE) > $@

assigned_reads_simlca.txt: kmers_rolling.txt $(TREE)
	$(TTIME) -o 4.2_read_assignment_simlca.log \
	$(ASSIGNMENT) -l -i $< -n $(TREE) > $@

5.1_contigs_stats.log: index.fa.fai
	../../bin/contig_statistics.py -k $(K) -f index.fa.fai > $@

_main_log.log: index.fa.$(K).bit.klcp 5.1_contigs_stats.log \
	kmers_rolling.txt kmers_restarted.txt kmers_rolling_skipping.txt kmers_restarted_skipping.txt \
	assigned_reads.txt assigned_reads_simlca.txt
	du -sh *.fa.* | grep -v "fa.amb" | grep -v "fa.fai" > 5.2_index_size.log
	echo > _main_log.log
	date >> _main_log.log
	pwd >> _main_log.log
	echo >> _main_log.log

	tail -n +1 [0-9]*.log >> _main_log.log

_main_log.md: _main_log.log 
	cat _main_log.log | ../../bin/reformat_log.py > _main_log.md

clean:
	rm -f index.fa index.fa.* Makefile.generated *.log
	rm -fr index/
