.PHONY: test


test: testref.fa.bwt
	../prophyle_index -a 2 testref.fa query.fa

testref.fa.bwt: testref.fa
	bwa index testref.fa
