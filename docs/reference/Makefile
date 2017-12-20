.PHONY: all help clean

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

all:
	./save_help.sh

help: ## Print help message
	    @echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm *.txt

