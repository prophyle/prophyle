.PHONY: all clean

include ../conf.mk

GENLIBS=bacteria viruses plasmids
FNS=$(addprefix seqs/,$(addsuffix /.complete.1, $(GENLIBS)))

all: $(FNS)
	echo $(FNS)
	echo $(GENLIBS)

seqs/%/.complete.1:
	mkdir -p seqs
	$(PROP) download -d seqs $*

clean:
	rm -fr seqs
