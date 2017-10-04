.PHONY: all clean

include ../conf.mk

FA1=test_L.fa
FA2=test_R.fa

K=7 10

DIFFS1 = $(addsuffix .txt, $(addprefix __diff_L., $(K)))
DIFFS2 = $(addsuffix .txt, $(addprefix __diff_R., $(K)))

all: $(DIFFS1) $(DIFFS2)
	@for f in $^; do \
		if [[ -s $$f ]]; then \
			echo "file $$f is not empty"; \
			exit 1; \
		fi; \
	done

__diff_L.%.txt: _join_L.%.txt _inp_L.%.txt
	diff -c $^ | tee $@

__diff_R.%.txt: _join_R.%.txt _inp_R.%.txt
	diff -c $^ | tee $@

_join_L.%.txt: _intersect.%.txt _out_L.%.txt
	sort -m $^ > $@

_join_R.%.txt: _intersect.%.txt _out_R.%.txt
	sort -m $^ > $@

_intersect.%.txt: _intersect.%.fa
	$(F2K) -m a -i $< -k $* > $@

_out_L.%.txt: _out_L.%.fa
	$(F2K) -m a -i $< -k $* > $@

_out_R.%.txt: _out_R.%.fa
	$(F2K) -m a -i $< -k $* > $@

_inp_L.%.txt:
	$(F2K) -m a -i $(FA1) -k $* > $@

_inp_R.%.txt:
	$(F2K) -m a -i $(FA2) -k $* > $@

_out_L.%.fa _out_R.%.fa _intersect.%.fa:
	$(ASM) -k $* -i $(FA1) -i $(FA2) -o _out_L.$*.fa -o _out_R.$*.fa -x _intersect.$*.fa

clean:
	rm -f _*

