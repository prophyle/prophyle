.PHONY: all clean compile_assembler compile_index

all: compile_assembler compile_index

clean:
	$(MAKE) -C assembler clean
	$(MAKE) -C kmer-index clean

compile_assembler:
	@cd assembler && cmake .
	$(MAKE) -C assembler

compile_index:
	$(MAKE) -C kmer-index
