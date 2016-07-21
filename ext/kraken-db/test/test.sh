#!/bin/bash

TMP=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
KMER_N=1000
TAXA_N=100
TEST_DB=${DIR}/${TMP}/test_db.kdb
SORTED_DB=${DIR}/${TMP}/sorted.kdb
SORTED_IDX=${DIR}/${TMP}/sorted.idx
TAXA_MAP=${DIR}/${TMP}/taxa_map.dmp
TAXA_LIST=${DIR}/${TMP}/taxa_list.txt

mkdir -p ${TMP}

${DIR}/build_test_db -k ${KMER_N} -t ${TAXA_N} -o ${TEST_DB}
${DIR}/taxid_sorter -d ${TEST_DB} -o ${SORTED_DB} -i ${SORTED_IDX} -m ${TAXA_MAP} -l ${TAXA_LIST}
