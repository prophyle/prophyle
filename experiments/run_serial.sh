#! /usr/bin/env bash

cd "$(dirname "$0")"

for d in 0*/; do
	make -C $d $@
done

