#! /usr/bin/env bash -u -f -e -o pipefail

# download 1000 reads (from http://www.ebi.ac.uk/ena/data/view/PRJEB2186)
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR015/ERR015599/ERR015599_1.fastq.gz | gzip -d --stdout | head -n 4000 > reads.fq

# classify reads using ProPhyle
prophyle classify sparc_index/ reads.fq > read_assignment.sam