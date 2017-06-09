#! /usr/bin/env bash -u -f -e -o pipefail

K=25

# download SPARC assemblies
(cd assemblies && ./download_assemblies.sh)

# create a tree compatible with ProPhyle
(cd tree && ./create_sparc_tree.py)

# build a ProPhyle index
prophyle index -k $K -g assemblies ./tree/sparc.nw sparc_index
