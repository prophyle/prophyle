.PHONY: all clean

include ../conf.mk

K=3

all:
	$(1STEP) -k $(K) $(FA) reads.fq > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(IND) query -k $(K) -b $(FA) reads.fq > _obtained_qual.txt
	diff _obtained_qual.txt exp_results_qual.txt

clean:
	rm -f _*.txt index.fa.*
