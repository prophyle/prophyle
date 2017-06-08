#! /usr/bin/env bash -f -e -o pipefail

prophyle=./prophyle/prophyle.py

$prophyle download bacteria
$prophyle index -k 5 -s 0.05 -T ~/prophyle/bacteria.nw test_index

