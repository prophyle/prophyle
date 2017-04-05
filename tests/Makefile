A = $(wildcard testA*/.)
B = $(wildcard testB*/.)
AB = $(A) $(B)

.PHONY : all clean big $(AB) A B

V=1

all: A

A: $(A)

B: $(B)

$(AB):
	@echo
	@echo "========================================="
	@echo "$@"
	@echo "========================================="
	@echo
	$(MAKE) -C $@

clean:
	@for a in $(AB); do \
		$(MAKE) -C $$a clean; \
	done;
