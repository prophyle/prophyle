.PHONY: all clean

include ../conf.mk

K=31
tree=tree.nw
match=index_match.txt
CPP_ASS=$(PROP_DIR)/prophyle_assignment/prophyle_assignment
PY_ASS=$(PROP_DIR)/prophyle_assignment.py

all:
	$(CPP_ASS) -f sam -m h1 -A $(tree) $(K) $(match) > _test_c.1.sam
	samtools view _test_c.1.sam | sort -k 3 | sort >_test_c.2.sam
	samtools view -H _test_c.1.sam | grep -v "^@PG"> _test_c.h.sam
	#
	$(PY_ASS) -f sam -m h1 -A $(tree) $(K) $(match) >_test_py.1.sam
	samtools view \
		-x h2 -x c2 -x hf -x cf -x ln -x ii -x is \
		_test_py.1.sam | sort -k 3 | sort >_test_py.2.sam
	samtools view -H _test_py.1.sam | grep -v "^@PG" > _test_py.h.sam
	#
	# test sequences
	diff _test_c.2.sam _test_py.2.sam > diff_sam.seqs.txt
	# test headers
	diff _test_c.h.sam _test_py.h.sam > diff_sam.headers.txt
	#
	#
	$(CPP_ASS) -f kraken -m c1 -D $(tree) $(K) $(match) | sort -k 3 | sort -k 2 >_test_c.kra
	$(PY_ASS) -f kraken -m c1 -D $(tree) $(K) $(match) | sort -k 3 | sort -k 2 >_test_py.kra
	diff _test_c.kra _test_py.kra > diff_kra.txt

clean:
	rm -f _test_* diff_*.txt
