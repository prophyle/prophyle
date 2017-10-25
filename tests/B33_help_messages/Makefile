.PHONY: all clean helpmsg

include ../conf.mk

SUBCOMMANDS = $(shell $(PROP) -h 2>&1 | grep -E "^ {4}[^ ]" | awk '{print $$1}' | xargs echo)

all: helpmsg

helpmsg:
	# tricks:
	#   - " " is used to test prophyle without subcommands
	#   - exit $$(($$? - 2 )) subtract 2 from the error code (command without any params should exit with 2)
	for x in " " $(SUBCOMMANDS); do \
		echo; \
		echo "================================"; \
		echo "prophyle $$x"; \
		echo "================================"; \
		echo; \
		$(PROP) $$x || (exit $$(($$? - 2 )) ); \
		echo; \
		$(PROP) $$x -h; \
		echo; \
	done

clean:
