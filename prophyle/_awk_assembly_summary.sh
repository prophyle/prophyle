#! /bin/bash

awk -F "\t" -v OFS="\t" '$12=="Complete Genome" && $11=="latest"\
  {print $1, $7, $20}' assembly_summary.txt >ftpselection.tsv
cut -f 3 ftpselection.tsv | sed 's/ftp:\/\//http:\/\//g' |\
  awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}\
  {ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' >ftpfilepaths.tsv
cut -f 1,2 ftpselection.tsv >acc2taxid.tsv
