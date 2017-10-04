.PHONY: all clean
.NOTPARALLEL: _match*

include ../conf.mk

K=14 16

DIFFS = $(addsuffix .txt, $(addprefix __diff., $(K)))

all: $(DIFFS)
	@for f in $^; do \
		if [[ -s "$$f" ]]; then \
			echo "file $$f is not empty"; \
			exit 1; \
		fi; \
	done

__diff.%.txt: _match.rolling_window.%.txt _match.restarted_search.%.txt
	diff -c $^ | tee $@

_match.rolling_window.%.txt:
	$(1STEP) -u -k $* $(FA) $(FQ) > $@

_match.restarted_search.%.txt:
	$(1STEP) -k $* $(FA) $(FQ) > $@

clean:
	rm -f _* index.fa.*
