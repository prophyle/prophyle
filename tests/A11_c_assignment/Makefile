.PHONY: all clean cpp py

include ../conf.mk

K=31
tree=tree.nw
match=index_match.txt
CPP_ASS=$(PROP_DIR)/prophyle_assignment/prophyle_assignment
PY_ASS=$(PROP_DIR)/prophyle_assignment.py

all: cpp py
	# test sequences
	diff _test.2.cpp.sam _test.2.py.sam | tee _diff_sam.seqs.txt | head -n 20
	# test headers
	diff _test.h.cpp.sam _test.h.py.sam | tee _diff_sam.headers.txt | head -n 20
	# test kraken-like output
	diff _test.kraken.cpp.txt _test.kraken.py.txt | tee _diff_kra.txt | head -n 20

py:
	$(PY_ASS) -f sam -m h1 -A $(tree) $(K) $(match) >_test.1.py.sam
	$(PY_ASS) -f kraken -m c1 $(tree) $(K) $(match) | sort > _test.kraken.py.txt
	samtools view \
		-x h2 -x c2 -x hf -x cf -x ln -x ii -x is\
		-x ti -x sn -x ra \
		_test.1.py.sam | sort >_test.2.py.sam
	samtools view -H _test.1.py.sam | grep -v "^@PG" > _test.h.py.sam

cpp:
	$(CPP_ASS) -f sam -m h1 -A $(tree) $(K) $(match) > _test.1.cpp.sam
	$(CPP_ASS) -f kraken -m c1 -D $(tree) $(K) $(match) | sort > _test.kraken.cpp.txt
	samtools view \
		-x ti -x sn -x ra \
		_test.1.cpp.sam | sort >_test.2.cpp.sam
	samtools view -H _test.1.cpp.sam | grep -v "^@PG"> _test.h.cpp.sam

clean:
	rm -f _*
