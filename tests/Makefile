A = $(wildcard testA*/.)
B = $(wildcard testB*/.)
AB = $(A) $(B)

PROP=../prophyle/prophyle.py
SUBCOMMANDS = $(shell $(PROP) -h 2>1 | grep -E "^ {4}[^ ]" | awk '{print $$1}' | xargs echo)
$(info $(SUBCOMMANDS))

.PHONY : all clean big $(AB) A B src

V=1

all: A helpmsg

A: $(A)

B: $(B)

$(AB): src
	@echo
	@echo "========================================="
	@echo "$@"
	@echo "========================================="
	@echo
	$(MAKE) -C $@

src:
	$(MAKE) -C ../prophyle

clean:
	@for a in $(AB); do \
		$(MAKE) -C $$a clean; \
	done;

helpmsg:
	for x in $(SUBCOMMANDS); do \
		$(PROP) $$x -h; \
	done
