.PHONY: all clean

BWA=../../bin/bwa
EXK=../../bin/prophyle-index

FA=index.fa
K=10
diffs = $(wildcard ./diff*.txt)

all: 
	$(BWA) index $(FA) $(FQ) 
	$(EXK) index -k $(K) $(FA) $(FQ)
	mv $(FA).sa $(FA).sa.separate
	mv $(FA).$(K).bit.klcp $(FA).$(K).bit.klcp.separate

	$(EXK) index -s -k $(K) $(FA) 

	cmp $(FA).sa $(FA).sa.separate > diff_sa.txt
	cmp $(FA).$(K).bit.klcp $(FA).$(K).bit.klcp.separate > diff_klcp.txt

	@for f in $(diffs); do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done


clean:
	rm -f index.fa.* *.txt

