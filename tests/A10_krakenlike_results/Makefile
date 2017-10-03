.PHONY: all clean krakenformat se_vs_pe

include ../conf.mk

K=2
tree=tree.nw
index=_idx_test
reads=reads.fq
ref=ref_res.txt

all: krakenformat se_vs_pe

krakenformat:
	$(PROP) index -k 2 -g . $(tree) $(index)
	$(PROP) classify $(index) $(reads) -f kraken | cut -f3,5 > _assigned_kmers.lossless.txt
	$(PROP) classify $(index) $(reads) -L -f kraken | cut -f3,5 > _assigned_kmers.tielca.txt
	$(PROP) classify $(index) $(reads) -M | cut -f3,5 > _observed.kraken.txt
	#diff expected.kraken.txt _observed.kraken.txt | tee _diff_ref.txt

se_vs_pe:
	# PAIRED END
	$(PROP) classify $(index) $(reads) | cut -f1-5 | sort -k 3 | sort >_test_single.txt
	echo
	$(PROP) classify $(index) $(reads) $(reads) | cut -f1-5 | sort -k 3 | sort >_test_paired.txt
	echo
	diff _test_single.txt _test_paired.txt > _diff_paired.txt
	# CHECK DIFF
	@for f in diff_*.txt ; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

clean:
	rm -rf _*
