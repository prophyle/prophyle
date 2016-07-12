SUBDIRS := $(wildcard */.)
#SUBDIRS_C := "Bacteria/.complete"
COMPLETES := Plasmids/.complete Viruses/.complete Bacteria/.complete

.PHONY: clean all #$(SUBDIRS)
.DEFAULT_GOAL := all

all: $(COMPLETES)

Bacteria/.complete:
	$(MAKE) -C Bacteria

Viruses/.complete:
	$(MAKE) -C Viruses

Plasmids/.complete:
	$(MAKE) -C Plasmids

clean:
	@for a in $$(ls); do \
		if [ -d $$a ]; then \
			$(MAKE) -C $$a clean; \
		fi; \
	done;
