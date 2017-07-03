#! /usr/bin/env bash -f -e -o pipefail

samtools view _test_bam1.bam  | ./prophyle/prophyle.py analyze ./_index_test - _test_analyze.basic
samtools view _test_bam1.bam  | ./prophyle/prophyle.py analyze -N ./_index_test - _test_analyze.ncbi

