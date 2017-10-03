.PHONY: all clean

include ../conf.mk

K=12

all:
	$(1STEP) -k $(K) -u -t 1 $(FA) $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(IND) query -k $(K) -u -t 2 $(FA) $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(IND) query -k $(K) -u -t 4 $(FA) $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

	$(IND) query -k $(K) -u -t 10 $(FA) $(FQ) > _obtained.txt
	diff _obtained.txt exp_results.txt

clean:
	rm -f _*.txt index.fa.*
