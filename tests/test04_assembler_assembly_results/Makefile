.PHONY: all clean

SHELL=/bin/bash -o pipefail

K=4
NORM=../../bin/_fa_norm.py
ASM=../../bin/assembler

.SECONDARY:

all:
	$(ASM) -k $(K) -i input.fa -i input.fa -x _output_obtained.fa -o /dev/null -o /dev/null

	$(NORM) -i expected_output.fa > _output_expected.norm.fa
	$(NORM) -i _output_obtained.fa > _output_obtained.norm.fa

	diff _output_expected.norm.fa _output_obtained.norm.fa | tee _diff.txt

clean:
	rm -f _*
