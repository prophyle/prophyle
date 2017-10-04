.PHONY: all clean krakenformat se_vs_pe

include ../conf.mk

K=2
tree=tree.nw
index=_index
reads=reads.fa
ref=ref_res.txt

all: krakenformat #se_vs_pe

krakenformat: _index.complete
	$(PROP) classify $(index) $(reads) -f kraken  > _assigned_kmers.lossless.txt
	$(PROP) classify $(index) $(reads) -L -f kraken  > _assigned_kmers.tielca.txt
	$(PROP) classify $(index) $(reads) -M  > _observed.kraken.txt
	diff expected.kraken.txt _observed.kraken.txt | tee _diff_ref.txt

se_vs_pe: _index.complete
	# PAIRED END
	$(PROP) classify $(index) $(reads) | cut -f1-5 | sort -k 3 | sort >_test_single.txt
	$(PROP) classify $(index) $(reads) $(reads) | cut -f1-5 | sort -k 3 | sort >_test_paired.txt
	diff _test_single.txt _test_paired.txt > _diff_paired.txt
	# CHECK DIFF
	@for f in diff_*.txt ; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

_index.complete:
	$(PROP) index -A -k 2 -g fastas $(tree) $(index)
	touch $@

clean:
	rm -rf _*
