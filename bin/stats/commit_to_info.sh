#! /usr/bin/env bash

#set -e
#set -o pipefail
#set -u
#set -f

readonly PROGNAME=$(basename $0)
readonly PROGDIR=$(dirname $0)
readonly -a ARGS=("$@")
readonly NARGS="$#"

case "$#" in
	1)
		commit="$1"
		;;
	*)
		s=`basename $0`
		echo "print info for a given ProPhyle commit"
		echo "usage: $s commit"
		exit 1
		;;
esac

echo "Processing commit $commit" >&2

(
git checkout "$commit"
) > /dev/null 2> /dev/null

commit=$(git rev-parse HEAD | head -c 6)
md5=$(find "$PROGDIR/../.." -type f -name '*.py' -o -name '*.nw' -o -name '*.c' -o -name '*.cpp' -o -name '*.h' -o -name 'Makefile' | sort | xargs cat | md5sum | cut -d ' ' -f 1)
version=$(find "$PROGDIR/../.." -type f -name '*.py' | xargs cat | grep -E "^VERSION\s+=" | perl -pe "s/.*=\s+\"(\d+(.\d+)+)\".*/\1/g")

(
git checkout master
) > /dev/null 2> /dev/null

printf "$commit\t$md5\t$version\n"

