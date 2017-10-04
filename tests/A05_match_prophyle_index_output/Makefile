.PHONY: all clean normal with_seqs

include ../conf.mk
FA=index.fa

K=3

all: normal with_seqs

normal: index.complete
	$(IND) query -k $(K) $(FA) reads.fq > _obtained.normal.txt
	diff -c expected.normal.txt _obtained.normal.txt

with_seqs: index.complete
	$(IND) query -k $(K) -b $(FA) reads.fq > _obtained.with_seqs.txt
	diff -c expected.with_seqs.txt _obtained.with_seqs.txt

index.complete:
	$(BWA) index $(FA)
	$(IND) build -k $(K) $(FA)
	touch $@

clean:
	rm -f _* index.fa.* *.complete
