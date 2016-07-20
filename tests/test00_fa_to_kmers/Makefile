.PHONY: all clean

F2K=../../bin/_fa_to_kmers.py -k 2 -i fasta.fa

all: gen_a.txt gen_c.txt gen_r.txt gen_f.txt
	diff exp_a.txt gen_a.txt
	diff exp_c.txt gen_c.txt
	diff exp_f.txt gen_f.txt
	diff exp_r.txt gen_r.txt

gen_%.txt:
	$(F2K) -m $* > $@

clean:
	rm -f gen*.txt