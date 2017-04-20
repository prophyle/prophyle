#! /usr/bin/env bash

set -e -f -o pipefail

cd $(dirname $0)

for x in main download index classify ; do
	fn="$x.txt"
	echo "$fn"

	if [[ "$x" == "main" ]] ; then
		x=""
	fi

	echo "$ prophyle $x -h" > $fn
	echo >> $fn
	prophyle $x -h >> $fn
done

# TODO: solve problem with dashes and output redirection
# TODO: change assignment.py to prophyle-assignment

# for x in assembler index assignment ; do
# 	fn="prophyle_$x.txt"
#		com="prophyle-$x"
# 	echo "$fn"
# 	echo "$ $com -h" > $fn
# 	echo >> $fn
# 	../../bin/$com >> $fn
# done
