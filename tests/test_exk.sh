#! /usr/bin/env bash

set -eu
set -o verbose
set -o pipefail

cd "$(dirname "$0")"

../bin/exk index -k 7 bwa_index/index.fa

../bin/exk match -k 7 -v -u bwa_index/index.fa ./reads/simulation_bacteria.1000.fq > match_rolling.txt
../bin/exk match -k 7 -v bwa_index/index.fa ./reads/simulation_bacteria.1000.fq > match_restarted.txt

diff match_rolling.txt match_restarted.txt

../bin/exk match -s -k 7 -v -u bwa_index/index.fa ./reads/simulation_bacteria.1000.fq > match_rolling_skip.txt
../bin/exk match -s -k 7 -v bwa_index/index.fa ./reads/simulation_bacteria.1000.fq > match_restarted_skip.txt

diff match_rolling_skip.txt match_restarted_skip.txt
