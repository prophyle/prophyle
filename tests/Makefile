A = $(wildcard testA*/.)
B = $(wildcard testB*/.)
AB = $(A) $(B)

.PHONY : all clean big $(AB) A B src

V=1

all: A

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
