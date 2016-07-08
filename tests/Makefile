SUBDIRS := $(wildcard test*/.)
V=1

.PHONY : all clean $(SUBDIRS)

all : $(SUBDIRS)

$(SUBDIRS): 
	$(MAKE) -C $@

clean:
	@for a in $(SUBDIRS); do \
		$(MAKE) -C $$a clean; \
	done;	