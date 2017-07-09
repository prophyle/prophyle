#! /usr/bin/env bash -u -f -e -o pipefail

K=25

# build a ProPhyle index
prophyle index -k $K -g assemblies -A ./tree/sparc.core_genes.nw sparc_index

