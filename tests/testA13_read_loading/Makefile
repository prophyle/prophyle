.PHONY: all clean
.NOTPARALLEL: match*

include ../conf.mk

K=16

all:
	$(1STEP) -k $(K) $(FA) $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(1STEP) -k $(K) -r 10 $(FA) $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(1STEP) -k $(K) -r 1000 $(FA) $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(1STEP) -k $(K) -r 100000 $(FA) $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

clean:
	rm -f index.fa.* _*.txt
