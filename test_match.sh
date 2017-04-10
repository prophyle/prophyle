#! /usr/bin/env bash -f -e -o pipefail

head -n 4 tests/simulation_bacteria.1000.fq | \
prophyle classify test_idx -

