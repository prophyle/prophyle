include ../conf.mk

K=12
T=1 2 3 4 5 6 7 8 9 10 20 50 100

.PHONY: all clean $(T)

all: $(addprefix t,$T)

t%: _index.%.complete
	$(IND) query -k $(K) -u -t $* _index.$*.fa $(FQ) > _obtained.$*.txt
	diff -c _obtained.$*.txt expected.txt

_index.%.complete:
	$(BWA) index -p _index.$*.fa $(FA)
	$(IND) build -k $(K) _index.$*.fa
	touch $@

clean:
	rm -f _*
