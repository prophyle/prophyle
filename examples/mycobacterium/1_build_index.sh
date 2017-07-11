#! /usr/bin/env bash -u -f -e -o pipefail

K=31

# download genomes
prophyle download bacteria

# build a ProPhyle index
prophyle index -k $K ~/prophyle/bacteria.nw@77643 mycobacterium_index
