#!/bin/bash

USAGE="kraken_assembler -t taxa_map_file"

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
	-t|--taxa-map)
	TAXA_MAP="$2"
	shift
	;;
	-d|--db)
	SORTED_DB="$2"
	shift
	;;
	-i|--idx)
	SORTED_IDX="$2"
	shift
	;;
	*)
	echo $USAGE
	exit(-1)
	;;
esac
shift # past argument or value
done
echo FILE EXTENSION  = "${EXTENSION}"
echo SEARCH PATH     = "${SEARCHPATH}"
echo LIBRARY PATH    = "${LIBPATH}"
echo "Number files in SEARCH PATH with EXTENSION:" $(ls -1 "${SEARCHPATH}"/*."${EXTENSION}" | wc -l)
if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

# while read line
# do
# 	../../bin/assembler \
# 		-i <() \
# 		-o /dev/null -x test/0.fa -k 31
# done < $1
