A = $(wildcard A*/.)
B = $(wildcard B*/.)
C = $(wildcard C*/.)

ABC = $(A) $(B) $(C)

$(info $(SUBCOMMANDS))

.PHONY: all clean big $(ABC) A B C compile parallel

.NOTPARALLEL: B00_compile

V=1

all: A B

parallel:
	$(MAKE) -j | tee __tests_parallel.log

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
	rm -f __*
	@for a in $(ABC); do \
		$(MAKE) -C $$a clean; \
	done;

