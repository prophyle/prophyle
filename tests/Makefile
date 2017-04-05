SUBDIRS := $(wildcard testA*/.)
V=1

.PHONY : all clean $(SUBDIRS)

all : $(SUBDIRS)

$(SUBDIRS):
	@echo
	@echo "========================================="
	@echo "$@"
	@echo "========================================="
	@echo
	$(MAKE) -C $@

clean:
	@for a in $(SUBDIRS); do \
		$(MAKE) -C $$a clean; \
	done;
