#! /usr/bin/env bash

/usr/bin/env pypy fai2taxid.py "$1" taxonomy/gi_taxid_nucl.dmp > 1_seq_taxid.txt
