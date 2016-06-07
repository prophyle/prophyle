SUBDIRS := $(wildcard */.)

.PHONY: clean all $(SUBDIRS)

all: .complete

.complete: $(SUBDIRS)
	bash -c "find . -name '*.fna' -o -name '*.ffn' | parallel --no-notice --verbose samtools faidx {}"
	touch .complete

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	@for a in $$(ls); do \
		if [ -d $$a ]; then \
			$(MAKE) -C $$a clean; \
		fi; \
	done;	