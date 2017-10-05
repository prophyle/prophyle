A = $(wildcard A*/.)
B = $(wildcard B*/.)
C = $(wildcard C*/.)

ABC = $(A) $(B) $(C)

$(info $(SUBCOMMANDS))

.PHONY: all clean big $(ABC) A B C compile

.NOTPARALLEL: B00_compile

V=1

all: A B

A: $(A)

B: $(B)

C: $(C)

$(ABC): compile
	@echo
	@echo "========================================="
	@echo "$@"
	@echo "========================================="
	@echo
	$(MAKE) -C $@

compile:
	$(MAKE) -C ../prophyle

clean:
	@for a in $(ABC); do \
		$(MAKE) -C $$a clean; \
	done;

