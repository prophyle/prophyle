#! /usr/bin/env bash

set -u
set -f
set -e
set -o pipefail

K=31

# download genomes
../../prophyle/prophyle.py download bacteria

# build a ProPhyle index
../../prophyle/prophyle.py index -k $K ~/prophyle/bacteria.nw@77643 mycobacterium_index
