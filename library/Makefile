SUBDIRS := $(wildcard */.)
#SUBDIRS_C := "Bacteria/.complete"
COMPLETES := Bacteria/.complete Plasmids/.complete Viruses/.complete HMP/.complete


.PHONY: clean all #$(SUBDIRS)
.DEFAULT_GOAL := all

all: $(COMPLETES)

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
