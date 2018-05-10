.PHONY: all clean index index_noprop

include ../conf.mk

K=5

all: index index_noprop

index_noprop:
	$(PROP) index -R -T -k $(K) -g lib tree.nw $@
	find $@/propagation -name "*.fa" -print | xargs -I {} -P 4 bash -c 'diff -c <($(F2K) -i "{}" -k $(K) -m c -f fa) expected_noprop/`basename "{}"`'

index:
	$(PROP) index -T -k $(K) -g lib tree.nw index
	find $@/propagation -name "*.fa" -print | xargs -I {} -P 4 bash -c 'diff -c <($(F2K) -i "{}" -k $(K) -m c -f fa) expected/`basename "{}"`'

clean:
	rm -rf index*
