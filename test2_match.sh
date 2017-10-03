#! /usr/bin/env bash

set -f
set -e
set -o pipefail

echo
echo " =========== 1 ============"
echo

set -o verbose
./prophyle/prophyle.py classify -P _index_test tests/simulation_bacteria.1000.fq | samtools view -b > _test_bam1.bam
set +o verbose

echo
echo " =========== 2 ============"
echo

set -o verbose
./prophyle/prophyle.py classify -P _index_test tests/simulation_bacteria.1000.fq tests/simulation_bacteria.1000.fq | samtools view -b > _test_bam2.bam
set +o verbose

