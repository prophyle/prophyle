#! /usr/bin/env bash

set -f
set -e
set -o pipefail

set -v

# =========== 1 ============

./prophyle/prophyle.py classify -P _index_test tests/simulation_bacteria.1000.fq | samtools view -b > _test_bam1.bam

# =========== 2 ============

./prophyle/prophyle.py classify -P _index_test tests/simulation_bacteria.1000.fq tests/simulation_bacteria.1000.fq | samtools view -b > _test_bam2.bam
