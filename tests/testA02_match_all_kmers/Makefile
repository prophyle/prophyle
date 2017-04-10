.PHONY: all clean
.NOTPARALLEL: prophyle_index_output*

include ../conf.mk

FA=index.fa
K=2 3 4 5 6 7 8

diffs = $(addsuffix .txt, $(addprefix diff., $(K)))

all: $(diffs)
	@for f in $^; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

clean:
	rm -f *.fq *.txt index.fa.*

kmers_all.%.fq:
	$(AK) -k $* -f fq > $@

kmers_all.%.txt:
	$(AK) -k $* -f txt > $@

prophyle_index_output.%.txt: kmers_all.%.fq
	$(1STEP) -v -k $* $(FA) $< > $@

kmers_matched.%.txt: prophyle_index_output.%.txt
	cat $< | (grep -B 1 '^[1-9]' || true) | (grep "#" || true) | sed 's/\#//' | sort > $@

kmers_expected.%.txt: $(FA)
	$(F2K) -i $< -k $* > $@

diff.%.txt: kmers_expected.%.txt kmers_matched.%.txt
	diff $^ | tee $@
