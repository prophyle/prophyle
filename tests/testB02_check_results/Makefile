.PHONY: all clean

include ../conf.mk

K=2
tree=tree.nw
index=test_index
reads=reads.fq
ref=ref_res.tsv
res=test_res.tsv
diff_res=diff_res.txt
diffs = $(wildcard ./diff*.txt)

all:
	prophyle index -k 2 -g . $(tree) $(index)
	prophyle classify $(index) $(reads) -L -f kraken | grep ^[C\t] | cut -f3,5 > $(res)
	cmp $(ref) $(res) > $(diff_res)
	@for f in $(diffs); do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

clean:
	rm -rf $(index) $(res) $(diff_res)
