#! /usr/bin/env bash

for x in main init index classify ; do
	fn="$x.txt"
	echo "$fn"

	if [[ "$x" == "main" ]] ; then
		x=""
	fi

	echo "$ prophyle $x -h" > $fn
	echo >> $fn
	../../bin/prophyle.py $x -h >> $fn
done
