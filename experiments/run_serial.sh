#! /usr/bin/env bash

for d in 0*/; do
	make -C $d $@
done

