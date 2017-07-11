#! /usr/bin/env bash

set -u
set -f
set -e
set -o pipefail

K=25

# build a ProPhyle index
../../prophyle/prophyle.py index -k $K -g assemblies -A ./tree/sparc.core_genes.nw sparc_index

