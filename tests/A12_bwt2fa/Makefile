include ../conf.mk

.PHONY: all clean $(T)

all: _index.complete _bwt2fa 
	diff -c $(FA) bwt2fa_index.fa

_bwt2fa: _index.complete
	$(IND) bwt2fa index.fa bwt2fa_index.fa

_index.complete:
	$(BWA) index $(FA)

clean:
	rm -f _*
