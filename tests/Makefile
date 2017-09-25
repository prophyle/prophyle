A = $(wildcard testA*/.)
B = $(wildcard testB*/.)
C = $(wildcard testB*/.)

ABC = $(A) $(B) $(C)

$(info $(SUBCOMMANDS))

.PHONY : all clean big $(ABC) A B C src

V=1

all: A B

A: $(A)

B: $(B)

C: $(C)

$(ABC): src
	@echo
	@echo "========================================="
	@echo "$@"
	@echo "========================================="
	@echo
	$(MAKE) -C $@

src:
	$(MAKE) -C ../prophyle

clean:
	@for a in $(ABC); do \
		$(MAKE) -C $$a clean; \
	done;

