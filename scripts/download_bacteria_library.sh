#!/bin/bash
cd "$(dirname "$0")"

wget http://monge.univ-mlv.fr/~brinda/library/library__full_bacteria.tar.gz
tar xvf library__full_bacteria.tar.gz
mv -f library_full_bacteria/Bacteria ../libraries
