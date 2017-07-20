#! /usr/bin/env bash -u -f -e -o pipefail

# download SPARC tree and assemblies
(cd tree && ./download.sh)
(cd assemblies && ./download.sh)
(cd reads && ./download.sh)
