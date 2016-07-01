.PHONY: all clean
.NOTPARALLEL: exk_output*

SHELL=/bin/bash -o pipefail

.SECONDARY:

M=../../bin/exk_match_1step.py
AK=../../bin/_all_kmers.py
FK=../../bin/_fa_to_kmers.py

FA=index.fa
K=2 3 4 5 6 7 8

targets = $(addprefix diff., $(K))

all: $(addsuffix .txt, $(targets))

clean:
	rm -f *.fq *.txt index.fa.*

kmers_all.%.fq:
	$(AK) -k $* -f fq > $@

kmers_all.%.txt:
	$(AK) -k $* -f txt > $@

exk_output.%.txt: kmers_all.%.fq
	$(M) -k $* $(FA) $< > $@

kmers_matched.%.txt: exk_output.%.txt
	cat $< | grep -B 1 '^[1-9]' | grep "#" | sed 's/\#//' | sort > $@

kmers_expected.%.txt: $(FA)
	$(FK) -i $< -k $* > $@

diff.%.txt: kmers_expected.%.txt kmers_matched.%.txt
	diff $^ | tee $@
