.PHONY: all clean helpmsg

include ../conf.mk

SUBCOMMANDS = $(shell $(PROP) -h 2>&1 | grep -E "^ {4}[^ ]" | awk '{print $$1}' | xargs echo)

all: helpmsg

helpmsg:
	# " " here is used to test prophyle without subcommands
	for x in " " $(SUBCOMMANDS); do \
		echo; \
		echo "================================"; \
		echo "prophyle $$x"; \
		echo "================================"; \
		echo; \
		$(PROP) $$x; \
		echo; \
		$(PROP) $$x -h; \
		echo; \
	done

clean:
