.PHONY: all clean
.NOTPARALLEL: match*

SHELL=/bin/bash -o pipefail

.SECONDARY:

M=../../bin/exk_match_1step.py

FA=index.fa
FQ=simulation_bacteria.1000.fq
K=14

targets = $(addsuffix .txt, $(K))

all: \
	$(addprefix diff., $(targets)) \
	$(addprefix diff_s., $(targets)) \

match_rol_s.%.txt:
	$(M) -k $* $(FA) $(FQ) > $@

match_rol.%.txt:
	$(M) -k $* $(FA) $(FQ) > $@

match_res_s.%.txt:
	$(M) -k $* $(FA) $(FQ) > $@

match_res.%.txt:
	$(M) -k $* $(FA) $(FQ) > $@

diff_s.%.txt: match_rol_s.%.txt match_res_s.%.txt
	diff $^ | tee $@

diff.%.txt: match_rol.%.txt match_res.%.txt
	diff $^ | tee $@

clean:
	rm -f index.fa.* *.txt
