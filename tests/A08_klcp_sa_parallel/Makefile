.PHONY: all clean

include ../conf.mk

K=10
diffs = $(wildcard ./diff*.txt)

all: 
	$(BWA) index $(FA) $(FQ) 
	$(IND) build -k $(K) $(FA) $(FQ)
	mv $(FA).sa $(FA).sa.separate
	mv $(FA).$(K).klcp $(FA).$(K).klcp.separate

	$(IND) build -s -k $(K) $(FA) 

	cmp $(FA).sa $(FA).sa.separate > diff_sa.txt
	cmp $(FA).$(K).klcp $(FA).$(K).klcp.separate > diff_klcp.txt

	@for f in $(diffs); do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done


clean:
	rm -f index.fa.* *.txt

