SUBDIRS := $(wildcard */.)
#SUBDIRS_C := "Bacteria/.complete"
COMPLETES := Plasmids/.complete Viruses/.complete Bacteria/.complete HMP/.complete

.PHONY: clean all #$(SUBDIRS)
.DEFAULT_GOAL := all

all: $(COMPLETES) .complete

.complete: $(COMPLETES)
	bash -c "find . -name '*.fna' -o -name '*.ffn' | parallel --no-notice --verbose samtools faidx {}"
	touch .complete

Bacteria/.complete:
	$(MAKE) -C Bacteria

Viruses/.complete:
	$(MAKE) -C Viruses

Plasmids/.complete:
	$(MAKE) -C Plasmids

HMP/.complete:
	$(MAKE) -C HMP

clean:
	@for a in $$(ls); do \
		if [ -d $$a ]; then \
			$(MAKE) -C $$a clean; \
		fi; \
	done;
