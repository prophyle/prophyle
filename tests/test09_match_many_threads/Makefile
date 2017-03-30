.PHONY: all clean

FQ=../simulation_bacteria.1000.fq

all:
	../../bin/1step_match.py -k 12 -u -t 1 index.fa $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	../../bin/prophyle-index query -k 12 -u -t 2 index.fa $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	../../bin/prophyle-index query -k 12 -u -t 4 index.fa $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	../../bin/prophyle-index query -k 12 -u -t 10 index.fa $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

clean:
	rm -f _*.txt
