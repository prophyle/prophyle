.PHONY: all clean

all:
	../../bin/1step_exk_match.py -k 3 ref.fa reads.fq > _obtained.txt
	diff _obtained.txt exp_results.txt

clean:
	rm -f _*.txt