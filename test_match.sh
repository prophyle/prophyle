#! /usr/bin/env bash -f -e -o pipefail

echo
echo " =========== 1 ============"
echo

./prophyle/prophyle.py classify -P _index_test tests/simulation_bacteria.1000.fq | samtools view -b > _test_bam1.bam

echo
echo " =========== 2 ============"
echo

./prophyle/prophyle.py classify -P _index_test tests/simulation_bacteria.1000.fq tests/simulation_bacteria.1000.fq | samtools view -b > _test_bam2.bam

