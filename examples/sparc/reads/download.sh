#! /usr/bin/env bash

set -f
set -o pipefail

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR015/ERR015599/ERR015599_1.fastq.gz | gzip -d --stdout | head -n 4000 > reads.fq

