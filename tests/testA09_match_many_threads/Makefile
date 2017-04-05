.PHONY: all clean

include ../conf.mk

FQ=../simulation_bacteria.1000.fq

all:
	$(1STEP) -k 12 -u -t 1 index.fa $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(IND) query -k 12 -u -t 2 index.fa $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(IND) query -k 12 -u -t 4 index.fa $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(IND) query -k 12 -u -t 10 index.fa $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

clean:
	rm -f _*.txt
