#! /usr/bin/env bash -u -f -e -o pipefail

K=25

(cd assemblies && ./download_assemblies.sh)

(cd tree && ./create_sparc_tree.py)

prophyle index -k $K -g assemblies ./tree/sparc.nw sparc_index

