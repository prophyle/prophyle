#! /usr/bin/env bash
set -u
set -f
set -e
set -o pipefail

# classify reads using ProPhyle
../../prophyle/prophyle.py classify sparc_index/ reads/reads.fq > read_assignment.sam
