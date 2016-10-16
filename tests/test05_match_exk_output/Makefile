.PHONY: all clean

all:
	../../bin/1step_match.py -k 3 ref.fa reads.fq > _obtained.txt
	diff _obtained.txt exp_results.txt

	../../bin/prophyle-index match -k 3 -b ref.fa reads.fq > _obtained_qual.txt
	diff _obtained_qual.txt exp_results_qual.txt

clean:
	rm -f _*.txt