#! /usr/bin/env bash
set -f
set -e
set -o pipefail

./prophyle/prophyle.py download bacteria
./prophyle/prophyle.py index -k 5 -s 0.05 -T ~/prophyle/bacteria.nw _index_test

