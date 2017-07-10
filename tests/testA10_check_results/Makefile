.PHONY: all clean

include ../conf.mk

K=2
tree=tree.nw
index=_idx_test
reads=reads.fq
ref=ref_res.txt

all:
	$(PROP) index -k 2 -g . $(tree) $(index)
	$(PROP) classify $(index) $(reads) -L -f kraken | grep ^[C\t] | cut -f3,5 >_test_res.txt
	diff $(ref) _test_res.txt > diff_ref.txt
	# PAIRED END
	$(PROP) classify $(index) $(reads) >_test_single.txt
	$(PROP) classify $(index) $(reads) $(reads) >_test_paired.txt
	diff _test_single.txt _test_paired.txt > diff_paired.txt
	# CHECK DIFF
	@for f in diff_*.txt ; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

clean:
	rm -rf $(index) _test_*.txt diff_*.txt
