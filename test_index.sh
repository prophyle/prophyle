#! /usr/bin/env bash -f -e -o pipefail -x

./prophyle/prophyle.pydownload bacteria
./prophyle/prophyle.py index -k 5 -s 0.05 -T ~/prophyle/bacteria.nw _test_index

