#! /usr/bin/env bash -f -e -o pipefail


for i in $(seq 2); do
	echo
	echo " =========== $i ============"
	echo

	samtools view _test_bam${i}.bam  | ./prophyle/prophyle.py analyze -N ./_index_test _test_analyze${i}.basic -
	samtools view _test_bam${i}.bam  | ./prophyle/prophyle.py analyze -N ./_index_test _test_analyze${i}.ncbi -

done
