.PHONY: all clean compile_assembler compile_index #compile_ext

all: compile_assembler compile_index #compile_ext

clean:
	$(MAKE) -C assembler clean
	$(MAKE) -C kmer-index clean
	$(MAKE) -C kraken-db clean

compile_assembler:
	$(MAKE) -C assembler

compile_index:
	$(MAKE) -C kmer-index

#compile_ext:
#	$(MAKE) -C ../ext
