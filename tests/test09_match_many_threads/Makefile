.PHONY: all clean

all:
	../../bin/1step_match.py -k 12 -u -t 1 index.fa simulation_bacteria.1000.fq > _obtained.txt
	diff _obtained.txt exp_results.txt

	../../bin/prophyle-index match -k 12 -u -t 2 index.fa simulation_bacteria.1000.fq > _obtained.txt
	diff _obtained.txt exp_results.txt

	../../bin/prophyle-index match -k 12 -u -t 4 index.fa simulation_bacteria.1000.fq > _obtained.txt
	diff _obtained.txt exp_results.txt

	../../bin/prophyle-index match -k 12 -u -t 10 index.fa simulation_bacteria.1000.fq > _obtained.txt
	diff _obtained.txt exp_results.txt

clean:
	rm -f _*.txt
