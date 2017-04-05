.PHONY: all clean

FA=index.fa
FQ=../simulation_bacteria.1000.fq
K=10

all:
	../../bin/1step_match.py -k $(K) $(FA) $(FQ) > _exk_output.txt

	# test if number of reads is correct
	$(eval L=$(shell cat $(FQ) | grep "^@" | wc -l))
	test `awk '{print $$2}' ./_exk_output.txt | sort | uniq | wc -l` -eq $(L)

	# test if reported contigs are correct
	cat _exk_output.txt | ./extract_reported_contigs.py | sort | uniq | grep -v "^0$$" > _reported.txt
	cat index.fa.fai | sed "s/_/	/" | awk '{print $$1}' | sort | uniq > _present.txt

	cat _present.txt _reported.txt | sort | uniq > _present_reported.txt
	diff _present.txt _present_reported.txt

clean:
	rm -f _*.txt
