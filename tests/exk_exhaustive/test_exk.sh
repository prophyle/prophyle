#! /usr/bin/env bash

set -eu
set -o verbose
set -o pipefail

cd "$(dirname "$0")"

./kmer_gen.py 4
../../bin/exk_match_1step.py -k 4 bwa_index/test_seq.fa ./test_kmers.fq > match.txt
cat match.txt | grep -B 1 '^[1-9]' > results.txt
./test_equal.py
