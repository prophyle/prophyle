.PHONY: all clean
.NOTPARALLEL: _prophyle_index_output*

include ../conf.mk

K=$(shell seq 2 8)
DIFFS = $(addsuffix .txt, $(addprefix __diff., $(K)))

all: $(DIFFS)
	@for f in $^; do \
		if [[ -s "$$f" ]]; then \
			echo "file $$f is not empty"; \
			exit 1; \
		fi; \
	done

__diff.%.txt: _kmers.expected.%.txt _kmers.observed.%.txt
	diff -c $^ | tee $@

_kmers.expected.%.txt: $(FA)
	$(F2K) -i $< -k $* > $@

_kmers.observed.%.txt: _prophyle_index_output.%.txt
	cat "$<" | (grep -B 1 '^[1-9]' || true) | (grep "#" || true) | sed 's/\#//' | sort > $@

_prophyle_index_output.%.txt: _kmers.all.%.fq
	$(1STEP) -v -k $* $(FA) $< > $@

_kmers.all.%.fq:
	$(AK) -k $* -f fq > $@

clean:
	rm -f _* index.fa.*

