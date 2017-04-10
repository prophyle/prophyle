.PHONY: all clean

include ../conf.mk

all:
	$(1STEP) -k 3 index.fa reads.fq > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(IND) query -k 3 -b index.fa reads.fq > _obtained_qual.txt
	diff _obtained_qual.txt exp_results_qual.txt

clean:
	rm -f _*.txt index.fa.*
