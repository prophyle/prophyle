.PHONY: all clean

all: 
	prophyle download -d seqs viruses
	prophyle index seqs/test_viruses.nw test_index

clean:
	rm -fr test_index seqs
