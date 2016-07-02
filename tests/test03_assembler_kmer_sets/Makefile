.PHONY: all clean

SHELL=/bin/bash -o pipefail

ASM=../../bin/assembler
F2K=../../bin/_fa_to_kmers.py -m c

FA1=test_L.fa
FA2=test_R.fa

K= 2 7 10

diffs1 = $(addsuffix .txt, $(addprefix diff_L.k, $(K)))
diffs2 = $(addsuffix .txt, $(addprefix diff_R.k, $(K)))


.SECONDARY:

all: $(diffs1) $(diffs2)
	@for f in $^; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

clean:
	rm -f _*.fa *.txt

diff_L.k%.txt: _join_L.k%.txt _inp_L.k%.txt
	diff $^ | tee $@

diff_R.k%.txt: _join_R.k%.txt _inp_R.k%.txt
	diff $^ | tee $@

_join_L.k%.txt: _intersect.k%.txt _out_L.k%.txt
	sort -m $^ > $@

_join_R.k%.txt: _intersect.k%.txt _out_R.k%.txt
	sort -m $^ > $@

_intersect.k%.txt: _intersect.k%.fa
	$(F2K) -i $< -k $* > $@

_out_L.k%.txt: _out_L.k%.fa
	$(F2K) -i $< -k $* > $@

_out_R.k%.txt: _out_R.k%.fa
	$(F2K) -i $< -k $* > $@

_inp_L.k%.txt: test_L.fa
	$(F2K) -i $< -k $* > $@

_inp_R.k%.txt: test_R.fa
	$(F2K) -i $< -k $* > $@

_out1.%.fa _out2.%.fa _intersect.%.fa:
	$(ASM) -k $* -i $(FA1) -i $(FA2) -o _out1.$*.fa -o _out2.$*.fa -x _intersect.$*.fa
