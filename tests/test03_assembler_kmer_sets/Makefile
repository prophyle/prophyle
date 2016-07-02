.PHONY: all clean

SHELL=/bin/bash -o pipefail

ASM=../../bin/assembler
F2K=../../bin/_fa_to_kmers.py -m c

FA1=test1.fa
FA2=test2.fa

K=6 8

diffs1 = $(addsuffix .txt, $(addprefix diff1., $(K)))
diffs2 = $(addsuffix .txt, $(addprefix diff2., $(K)))


.SECONDARY:

all: $(diffs1) $(diffs2)
	@for f in $^; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

clean:
	rm _*.fa *.txt

diff1.%.txt: _join1.%.txt _inp1.%.txt
	diff $^ | tee $@

diff2.%.txt: _join2.%.txt _inp2.%.txt
	diff $^ | tee $@

_join1.%.txt: _intersect.%.txt _out1.%.txt
	sort -m $^ > $@

_join2.%.txt: _intersect.%.txt _out2.%.txt
	sort -m $^ > $@

_intersect.%.txt: _intersect.%.fa
	$(F2K) -i $< -k $* > $@

_out1.%.txt: _out1.%.fa
	$(F2K) -i $< -k $* > $@

_out2.%.txt: _out2.%.fa
	$(F2K) -i $< -k $* > $@

_inp1.%.txt: test1.fa
	$(F2K) -i $< -k $* > $@

_inp2.%.txt: test2.fa
	$(F2K) -i $< -k $* > $@

_out1.%.fa _out2.%.fa _intersect.%.fa:
	$(ASM) -k $* -i $(FA1) -i $(FA2) -o _out1.$*.fa -o _out2.$*.fa -x _intersect.$*.fa
