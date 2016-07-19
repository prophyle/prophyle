#!/bin/bash

USAGE="Usage: kraken_assembler -t taxa-map-file -d sorted-db \
		-i sorted-idx -f indexing-dir -k kmer-len"

while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-t|--taxa-map)
			TAXA_MAP=$2
			shift
		;;
		-d|--db)
			DB=$2
			shift
		;;
		-i|--idx)
			IDX=$2
			shift
		;;
		-f|--indexing-dir)
			INDEXING_DIR=$2
			if [[ "$INDEXING_DIR" != */ ]]; then
				INDEXING_DIR="$(INDEXING_DIR)/"
			fi
			shift
		;;
		-k)
			K=$2
			shift
		;;
		*)
			echo $USAGE
			exit 1
		;;
	esac
done



while read taxid; do
	../../bin/assembler \
		-i <(./get_kmers_by_taxid -t $(taxid) \
		-d $(DB) -i $(IDX) -m $(TAXA_MAP) ) \
		-o /dev/null -x $(INDEXING_DIR)$(taxid).fa -k $(K)
done < <(./taxa_list -m $(TAXA_MAP))
