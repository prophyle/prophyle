.PHONY: all renderfasta

include *.mk

#ASSEMBLER=../../bin/assembler
EXK=../../bin/exk
BUILD_FA=../../src/build_index.py
BWA=../../bin/bwa
SAMTOOLS?=samtools
NEWICK2MAKEFILE=../../bin/newick2makefile.py
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
TIME=../../bin/time
TTIME:=$(TIME)

all: index.fa.sa index.fa.$(K).bit.klcp _main_log.log _main_log.md

index/.complete: $(TREE)
	mkdir -p index

	$(NEWICK2MAKEFILE) \
	-n $(TREE) \
	-o ./index \
	-l ../../ \
	-k $(K) \
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
	$(TTIME) -o 2.2_bwa_pac2bwt.log \
	$(BWA) pac2bwt -d index.fa.pac index.fa.bwt

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
		-k $(K) -u -v index.fa $(READS) > $@

kmers_restarted.txt: $(READS) index.fa.sa \
	kmers_rolling.txt
	$(TTIME) -o 3.2a_matching_restarted.log \
	$(EXK) match -l 3.2b_matching_restarted.log \
		-k $(K) -v index.fa $(READS) > $@

kmers_rolling_skipping.txt: $(READS) index.fa.sa $(KLCP) \
	kmers_restarted.txt
	$(TTIME) -o 3.3a_matching_rolling_skipping.log \
	$(EXK) match -l 3.3b_matching_rolling_skipping.log \
		-k $(K) -u -v -s index.fa $(READS) > $@

kmers_restarted_skipping.txt: $(READS) index.fa.sa kmers_rolling.txt \
	kmers_rolling_skipping.txt
	$(TTIME) -o 3.4a_matching_restarted_skipping.log \
	$(EXK) match -l 3.4b_matching_restarted_skipping.log \
		-k $(K) -v -s index.fa $(READS) > $@

4.1_contigs_stats.log: index.fa.fai
	../../bin/contig_statistics.py -k $(K) -f index.fa.fai > 4.1_contigs_stats.log

_main_log.log: index.fa.$(K).bit.klcp 4.1_contigs_stats.log \
	kmers_rolling.txt kmers_restarted.txt kmers_rolling_skipping.txt kmers_restarted_skipping.txt
	du -sh *.fa.* | grep -v "fa.amb" | grep -v "fa.fai" > 4.2_index_size.log
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
