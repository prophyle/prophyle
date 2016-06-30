.PHONY: all clean

SHELL=/bin/bash -o pipefail

.SECONDARY:

AK=../../bin/_all_kmers.py
FK=../../bin/_fa_to_kmers.py

TESTED_K=2 3 4 5 6 7 8
targets = $(addprefix diff., $(TESTED_K))

all: $(addsuffix .txt, $(targets))

clean:
	rm -f *.fq *.txt *.fa

kmers_all.%.fq:
	$(AK) -k $* -f fq > $@

kmers_all.%.txt:
	$(AK) -k $* -f txt > $@

exk_output.%.txt: kmers_all.%.fq
	../../bin/exk_match_1step.py -k $* bwa_index/test_seq.fa $< > $@

kmers_matched.%.txt: exk_output.%.txt
	cat $< | grep -B 1 '^[1-9]' | grep "#" | sed 's/\#//' | sort > $@

kmers_expected.%.txt: bwa_index/test_seq.fa
	$(FK) -i $< -k $* > $@

diff.%.txt: kmers_expected.%.txt kmers_matched.%.txt
	diff $^ | tee $@
