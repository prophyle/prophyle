#!/bin/bash
cd "$(dirname "$0")"

wget http://monge.univ-mlv.fr/~brinda/library/library__full_bacteria.tar.gz
tar xvf library__full_bacteria.tar.gz
mkdir -p ../library
mv -f library_full_bacteria/Bacteria ../library
