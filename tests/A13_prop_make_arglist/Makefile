.PHONY: all clean index index_noprop lib

include ../conf.mk

K=5

all: index index_noprop

index_noprop: lib
	$(PROP) index -R -T -k $(K) -g lib tree.nw $@
	find $@/propagation -name '*.fa' | xargs -I {} bash -c 'diff -c <($(F2K) -i "{}" -k $(K) -m c -f fa) expected_noprop/`basename "{}"`'

index: lib
	$(PROP) index -T -k $(K) -g lib tree.nw $@
	find $@/propagation -name '*.fa' | xargs -I {} bash -c 'diff -c <($(F2K) -i "{}" -k $(K) -m c -f fa) expected/`basename "{}"`'

lib:
	cd lib && for ((i=6; i<=10000; i++)); do ln -s 5.fa "$$i.fa" ; done

clean:
	rm -rf index*
	for ((i=6; i<=10000; i++)); do rm -f "lib/$$i.fa" ; done
