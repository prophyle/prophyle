.PHONY: all clean

BWA=../../bin/bwa
PROPHYLE_INDEX=../../bin/prophyle-index

FA=index.fa
K=10
diffs = $(wildcard ./diff*.txt)

all: 
	$(BWA) index $(FA) $(FQ) 
	$(PROPHYLE_INDEX) build -k $(K) $(FA) $(FQ)
	mv $(FA).sa $(FA).sa.separate
	mv $(FA).$(K).klcp $(FA).$(K).klcp.separate

	$(PROPHYLE_INDEX) build -s -k $(K) $(FA) 

	cmp $(FA).sa $(FA).sa.separate > diff_sa.txt
	cmp $(FA).$(K).klcp $(FA).$(K).klcp.separate > diff_klcp.txt

	@for f in $(diffs); do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done


clean:
	rm -f index.fa.* *.txt

