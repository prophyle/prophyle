#! /usr/bin/env bash

readonly PROGDIR=$(dirname $0)
readonly ARGS="$@"

${PROGDIR}/prophyle/prophyle.py $ARGS
