.PHONY: all renderfasta

include *.mk

#ASSEMBLER=../../bin/assembler
EXK=../../bin/exk
BUILD_FA=../../src/build_index.py
BWA=bwa
NEWICK2MAKEFILE=../../bin/newick2makefile.py
FINAL_FA=../../bin/create_final_fasta.py

all: index.fa.bwt

index/:
	mkdir index

Makefile.generated:
	$(NEWICK2MAKEFILE) \
	-n $(TREE) \
	-o ./index \
	-l ../../ \
	-k $(K) \
	> Makefile.generated

index.fa: index/ Makefile.generated
	$(MAKE) -f Makefile.generated

	$(FINAL_FA) index > index.fa
	#@cat index/*.reduced.fa > index.fa

#index.fa:
#	$(BUILD_FA) \
#	-n $(TREE) \
#	-o ./index \
#	-l ../../ \
#	-k $(K) \
#
#	@cat index/*.reduced.fa > index.fa
#
#renderfasta:
#	@cat index/*.reduced.fa > index.fa

%.sa %.pac %.bwt %.amb %.ann: %
	$(BWA) index -a is $<

clean:
	rm -fr index index.fa index.fa.* Makefile.generated