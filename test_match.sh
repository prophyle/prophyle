#! /usr/bin/env bash -f -e -o pipefail

head -n 4 tests/simulation_bacteria.1000.fq | ./prophyle/prophyle.py classify _index_test -

