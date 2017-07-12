.PHONY: all clean

include ../conf.mk

K=31
tree=tree.nw
match=index_match.txt
CPP_ASS=$(PROP_DIR)/prophyle_assignment/prophyle_assignment
PY_ASS=$(PROP_DIR)/prophyle_assignment.py

all:
	$(CPP_ASS) -f sam -m h1 -A $tree $K $match | sort -k 3 | sort >_test_c.sam
	$(PY_ASS) -f sam -m h1 -A $tree $K $match | sort -k 3 | sort >_test_py.sam
	diff _test_c.sam _test_py.sam > diff_sam
	$(C_ASS) -k kraken -m c1 -D $tree $K $match | sort -k 3 | sort -k 2 >_test_c.kra
	$(PY_ASS) -k kraken -m c1 -D $tree $K $match | sort -k 3 | sort -k 2 >_test_py.kra
	diff _test_c.kra _test_py.kra > diff_kra
	# CHECK DIFF
	@for f in diff_*.txt ; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

clean:
	rm -f _test_* diff_*.txt
