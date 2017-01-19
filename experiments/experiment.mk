.PHONY: all renderfasta

V=1
SHELL:=/bin/bash -o pipefail

include *.mk

#ASSEMBLER=../../bin/assembler
EXK=../../bin/prophyle-index
BUILD_FA=../../src/build_index.py
BWA=../../bin/bwa
SAMTOOLS?=samtools
NEWICK2MAKEFILE=../../bin/newick2makefile.py
ASSIGNMENT=../../bin/assignment.py

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
	REP_PARAM:=-r
else
	REP_PARAM:=
endif

ifdef NONDEL
	FINAL_FA:=../../bin/create_final_fasta.py --nondel
else
	FINAL_FA:=../../bin/create_final_fasta.py
endif

ifdef NONPROP
	FINAL_FA:=../../bin/create_final_fasta.py --nondel
endif

all: index.fa.$(K).bit.klcp _main_log.log _main_log.md \
	assigned_reads.bam assigned_reads_simlca.bam

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
	$(MAKE) -f Makefile.generated V=1

	touch index/.complete

index.fa: index/.complete
	$(TTIME) -o 1.2_merging_fasta.log \
	$(FINAL_FA) index > $@

	# todo: add this to the main prophyle cli script
	touch $@.kmers.tsv
	echo "#file	no_kmers" >> $@.kmers.tsv
	cat index/*.count.tsv | grep -v "^#" | sort | uniq >> $@.kmers.tsv

index.fa.pac: index.fa
	$(TTIME) -o 2.1_bwa_fa2pac.log \
	$(BWA) fa2pac index.fa index.fa

index.fa.bwt: index.fa.pac
	$(TTIME) -o 2.2_bwa_pac2bwtgen.log \
	$(BWA) pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

	$(TTIME) -o 2.3_bwa_bwtupdate.log \
	$(BWA) bwtupdate index.fa.bwt

#index.fa.sa: index.fa.bwt
#	$(TTIME) -o 2.4_bwa_bwt2sa.log \
#	$(BWA) bwt2sa index.fa.bwt index.fa.sa

$(KLCP): index.fa index.fa.bwt
	$(TTIME) -o 2.5_klcp_sa.log \
	$(EXK) index -s -k $(K) $<

%.fai: %
	$(SAMTOOLS) faidx $<

kmers_rolling.txt: $(KLCP)
	$(TTIME) -o 3.1a_matching_rolling.log \
	$(EXK) match -b -l 3.1b_matching_rolling.log  \
		-k $(K) -u index.fa $(READS) > $@

kmers_restarted.txt: $(READS) $(KLCP) \
	kmers_rolling.txt
	$(TTIME) -o 3.2a_matching_restarted.log \
	$(EXK) match -b -l 3.2b_matching_restarted.log \
		-k $(K) index.fa $(READS) > $@

assigned_reads.bam: kmers_rolling.txt $(TREE)
	$(TTIME) -o 4.1_read_assignment.log \
	$(ASSIGNMENT) -i $< -n $(TREE) -k $(K) -f sam -a | $(SAMTOOLS) view -b > $@

assigned_reads_simlca.bam: kmers_rolling.txt $(TREE)
	$(TTIME) -o 4.2_read_assignment_simlca.log \
	$(ASSIGNMENT) -l -i $< -n $(TREE) -k $(K) -f sam -a -t | $(SAMTOOLS) view -b > $@

5.1_contigs_stats.log: index.fa.fai
	../../bin/contig_statistics.py -k $(K) -f index.fa.fai > $@

_main_log.log: index.fa.$(K).bit.klcp 5.1_contigs_stats.log \
	kmers_rolling.txt kmers_restarted.txt \
	assigned_reads.bam assigned_reads_simlca.bam
	du -sh *.fa.* | grep -v "fa.amb" | grep -v "fa.fai" > 5.2_index_size.log
	echo > _main_log.log
	date >> _main_log.log
	pwd >> _main_log.log
	echo >> _main_log.log

	tail -n +1 [0-9]*.log >> _main_log.log

_main_log.md: _main_log.log
	cat _main_log.log | ../../bin/reformat_log.py > _main_log.md

clean:
	rm -f index.fa index.fa.* Makefile.generated *.log *.bam
	rm -fr index/
