#! /usr/bin/env bash

# prophyle (main program)

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

# prophyle_index

fn="prophyle_index"
com="../../prophyle/prophyle_index/prophyle_index"

echo "$fn.txt"
echo "$ `basename $com` -h" > $fn.txt
echo >> $fn.txt
$com 2>> $fn.txt

for sub_com in build query ; do
	sub_fn="${fn}_${sub_com}.txt"
	echo $sub_fn
	echo "$ `basename $com` $sub_com -h" > $sub_fn
	echo >> $sub_fn
	$com $sub_com 2>> $sub_fn
done

# prophyle_assembler

fn="prophyle_assembler.txt"
com="../../prophyle/prophyle_assembler/prophyle_assembler"

echo "$fn"
echo "$ `basename $com` -h" > $fn
echo >> $fn
$com -h 2>> $fn

# prophyle_assignment

fn="prophyle_assignment.txt"
com="../../prophyle/prophyle_assignment.py"

echo "$fn"
echo "$ `basename $com` -h" > $fn
echo >> $fn
$com -h >> $fn
