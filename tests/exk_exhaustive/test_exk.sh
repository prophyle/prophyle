#! /usr/bin/env bash

set -eu
set -o verbose
set -o pipefail

cd "$(dirname "$0")"

../../bin/exk index -k 4  bwa_index/test_seq.fa

../../bin/exk match -k 4 -v -u bwa_index/test_seq.fa ./test_kmers.fq > match_rolling.txt
../../bin/exk match -k 4 -v bwa_index/test_seq.fa ./test_kmers.fq > match_restarted.txt

diff match_rolling.txt match_restarted.txt

../../bin/exk match -s -k 4 -v -u bwa_index/test_seq.fa ./test_kmers.fq > match_rolling_skip.txt
../../bin/exk match -s -k 4 -v bwa_index/test_seq.fa ./test_kmers.fq > match_restarted_skip.txt

diff match_rolling_skip.txt match_restarted_skip.txt

# cat match_restarted.txt | grep -B 1 '^[1-9]'