#! /usr/bin/env bash

set -f
set -o pipefail

#
# Download tree (see http://datadryad.org/resource/doi:10.5061/dryad.t55gq/1)
#
curl http://datadryad.org/bitstream/handle/10255/dryad.83423/SPARC.core_genes.tree?sequence=1 > sparc.core_genes.nw

