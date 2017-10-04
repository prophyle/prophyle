.PHONY: all clean krakenformat se_vs_pe LL L X K

include ../conf.mk

K=2
tree=tree.nw
index=_index
digrams=digrams.fa
ref=ref_res.txt

all: krakenformat

krakenformat: LL L X K #se_vs_pe
	for x in lossless assignment_lca kmer_lca krakenlike; do \
		diff -c expected.digrams.$$x.txt _assigned.digrams.$$x.txt \
			| tee _diff.digrams.$$x.txt \
			| head -n 20; \
	done

LL: _index.complete
	$(PROP) classify $(index) $(digrams) -f kraken | sort > _assigned.digrams.lossless.txt

L: _index.complete
	$(PROP) classify $(index) $(digrams) -f kraken -L | sort > _assigned.digrams.assignment_lca.txt

X: _index.complete
	$(PROP) classify $(index) $(digrams) -f kraken -X | sort > _assigned.digrams.kmer_lca.txt

K: _index.complete
	$(PROP) classify $(index) $(digrams) -M | sort > _assigned.digrams.krakenlike.txt

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
