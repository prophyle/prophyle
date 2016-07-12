.PHONY: all clean compile_assembler compile_index compile_ext compile_kraken

all: compile_assembler compile_index compile_ext compile_kraken

clean:
	$(MAKE) -C assembler clean
	$(MAKE) -C kmer-index clean
	$(MAKE) -C kraken-db clean

compile_assembler:
	@cd assembler && cmake .
	$(MAKE) -C assembler

compile_index:
	$(MAKE) -C kmer-index

compile_ext:
	$(MAKE) -C ../ext

compile_kraken:
	$(MAKE) -C kraken-db
