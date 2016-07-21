TMP=$1
KMER_N=1000
TAXA_N=100
TEST_DB=${TMP}/test_db.kdb
SORTED_DB=${TMP}/sorted.kdb
SORTED_IDX=${TMP}/sorted.idx
TAXA_MAP=${TMP}/taxa_map.dmp
TAXA_LIST=${TMP}/taxa_list.txt

mkdir -p $(TMP)

./build_test_db -k ${KMER_N} -t ${TAXA_N} -o ${TEST_DB}
./taxid_sorter -d ${TEST_DB} -o ${SORTED_DB} -i ${SORTED_IDX} -m ${TAXA_MAP} -l ${TAXA_LIST}

