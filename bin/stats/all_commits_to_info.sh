#! /usr/bin/env bash

git log --all --oneline | cut -f1 -d ' ' | xargs -n 1 ./commit_to_info.sh

