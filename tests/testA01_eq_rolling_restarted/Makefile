.PHONY: all clean
.NOTPARALLEL: match*

include ../conf.mk

K=14 16

#targets = $(addsuffix .txt, $(K))

diffs = $(addsuffix .txt, $(addprefix diff., $(K)))
diffs_s = $(addsuffix .txt, $(addprefix diff_s., $(K)))

all: $(diffs)
	@for f in $^; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

match_rol.%.txt:
	$(1STEP) -u -k $* $(FA) $(FQ) > $@

match_res.%.txt:
	$(1STEP) -k $* $(FA) $(FQ) > $@

diff.%.txt: match_rol.%.txt match_res.%.txt
	diff $^ | tee $@

clean:
	rm -f index.fa.* *.txt
