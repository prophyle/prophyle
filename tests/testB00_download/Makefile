.PHONY: all clean

GENLIBS=bacteria viruses plasmids
FNS=$(addprefix seqs/,$(addsuffix /.complete.1, $(GENLIBS)))

all: $(FNS)
	echo $(FNS)
	echo $(GENLIBS)

seqs/%/.complete.1:
	mkdir -p seqs
	prophyle download -d seqs $*

clean:
	rm -fr seqs
