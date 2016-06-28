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

#UNAME_S := $(shell uname -s)
#
#ifeq ($(UNAME_S),Darwin)
#	TIME?=gtime
#else
#	TIME?=time
#endif
#TTIME:=$(TIME) -v
#TTIME:=$(TIME)
TTIME:=time

all: index.fa.sa index.fa.$(K).bit.klcp _main_log.log _main_log.md

index/.complete: $(NEWICK2MAKEFILE) $(TREE)
	mkdir -p index

	$(NEWICK2MAKEFILE) \
	-n $(TREE) \
	-o ./index \
	-l ../../ \
	-k $(K) \
	> Makefile.generated

	$(TTIME) -o 1_kmer_propagation.log \
	$(MAKE) -f Makefile.generated

	touch index/.complete

index.fa: index/.complete
	$(TTIME) -o 2_merging_fasta.log \
	$(FINAL_FA) index > index.fa

7_contigs_stats.log: index.fa.fai
	../../bin/contig_statistics.py -k $(K) -f index.fa.fai > 7_contigs_stats.log

index.fa.pac: index.fa $(BWA)
	$(TTIME) -o 3.1_bwa_fa2pac.log \
	$(BWA) fa2pac index.fa index.fa

index.fa.bwt: index.fa.pac $(BWA)
	$(TTIME) -o 3.2_bwa_pac2bwt.log \
	$(BWA) pac2bwt index.fa.pac index.fa.bwt

	$(TTIME) -o 3.3_bwa_bwtupdate.log \
	$(BWA) bwtupdate index.fa.bwt

index.fa.sa: index.fa.bwt $(BWA)
	$(TTIME) -o 3.3_bwa_bwt2sa.log \
	$(BWA) bwt2sa index.fa.bwt index.fa.sa

index.fa.$(K).bit.klcp: index.fa index.fa.bwt index.fa.sa $(EXK)
	$(TTIME) -o 4_klcp.log \
	$(EXK) index -k $(K) $<

%.fai: %
	$(SAMTOOLS) faidx $<

kmers_rolling.txt: $(READS) index.fa.sa index.fa.$(K).bit.klcp $(EXK)
	$(TTIME) -o 5_matching_rolling.log \
	$(EXK) match -k $(K) -u -v index.fa $(READS) > kmers_rolling.txt

kmers_restarted.txt: $(READS) index.fa.sa kmers_rolling.txt $(EXK)
	$(TTIME) -o 6_matching_restarted.log \
	$(EXK) match -k $(K) -v index.fa $(READS) > kmers_restarted.txt

_main_log.log: index.fa.$(K).bit.klcp kmers_rolling.txt kmers_restarted.txt 7_contigs_stats.log
	du -sh *.fa.* | grep -v "fa.amb" > 8_index_size.log
	echo > _main_log.log
	date >> _main_log.log
	pwd >> _main_log.log
	echo >> _main_log.log

	tail -n +1 [0-9]*.log >> _main_log.log

_main_log.md: _main_log.log ../../bin/reformat_log.py
	cat _main_log.log | ../../bin/reformat_log.py > _main_log.md

clean:
	rm -f index.fa index.fa.* Makefile.generated *.log
	rm -fr index/
