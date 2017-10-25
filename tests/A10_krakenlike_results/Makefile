.PHONY: all clean krakenformat se_vs_pe LL L X K

include ../conf.mk

K=2
tree=data/tree.nw
index=_index
digrams=digrams.fa
ref=ref_res.txt

PROPC=$(PROP) classify -c conf.json

all: krakenformat

krakenformat: LL L X K #se_vs_pe
	for x in lossless assignment_lca kmer_lca krakenlike; do \
		diff -c expected.digrams.$$x.txt _assigned.digrams.$$x.txt \
			| tee __diff.digrams.$$x.txt \
			| head -n 20; \
	done

LL: _index.complete
	$(PROPC) -f kraken    $(index) $(digrams) | sort > _assigned.digrams.lossless.txt

L: _index.complete
	$(PROPC) -f kraken -L $(index) $(digrams) | sort > _assigned.digrams.assignment_lca.txt

X: _index.complete
	$(PROPC) -f kraken -X $(index) $(digrams) | sort > _assigned.digrams.kmer_lca.txt

K: _index.complete
	$(PROPC) -M           $(index) $(digrams) | sort > _assigned.digrams.krakenlike.txt

# TODO: finish this part of the test
se_vs_pe: _index.complete
	# PAIRED END
	$(PROPC) $(index) $(reads) | cut -f1-5 | sort -k 3 | sort >_test_single.txt
	$(PROPC) $(index) $(reads) $(reads) | cut -f1-5 | sort -k 3 | sort >_test_paired.txt
	diff -c _test_single.txt _test_paired.txt > __diff_paired.txt
	# CHECK DIFF
	@for f in diff_*.txt ; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

_index.complete:
	$(PROP) index -A -k 2 $(tree) $(index)
	touch $@

clean:
	rm -rf _*
