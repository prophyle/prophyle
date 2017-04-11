.PHONY: all clean

include ../conf.mk

K=4

all:
	$(ASM) -k $(K) -i input.fa -o _output_obtained.fa

	$(NORM) -i expected_output.fa > _output_expected.norm.fa
	$(NORM) -i _output_obtained.fa > _output_obtained.norm.fa

	diff _output_expected.norm.fa _output_obtained.norm.fa | tee _diff.txt

clean:
	rm -f _*
