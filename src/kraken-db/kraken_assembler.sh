#!/bin/bash

usage() {
	echo "Usage: $0 -l taxa-list -m taxa-map -d sorted-db -i sorted-idx -f indexing-dir -k kmer-len"
	exit 1
}


while getopts "l:m:d:i:f:k:" opt; do
	case $opt in
		l)
			TAXA_LIST=${OPTARG}
			;;
		m)
			TAXA_MAP=${OPTARG}
			;;
		d)
			DB=${OPTARG}
			;;
		i)
			IDX=${OPTARG}
			;;
		f)
			INDEXING_DIR=${OPTARG%/}
			;;
		k)
			K=${OPTARG}
			;;
		*)
			usage
			;;
	esac
done

if [ -z "${TAXA_LIST}" ] || [ -z "${TAXA_MAP}" ] || [ -z "${DB}" ] || \
	[ -z "${IDX}" ] || [ -z "${INDEXING_DIR}" ] || [ -z "${K}" ]
then
	usage
fi

mkdir -p ${INDEXING_DIR}
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while read taxid; do
	${DIR}/assembler \
		-i <(${DIR}/kraken/get_kmers_by_taxid -t ${taxid} \
		-d ${DB} -i ${IDX} -m ${TAXA_MAP}) \
		-o /dev/null -x ${INDEXING_DIR}/${taxid}.fa -k ${K}
done < ${TAXA_LIST}
