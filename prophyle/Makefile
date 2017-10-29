.PHONY: all clean compile_assembler compile_index submodules

# BWA submodule version
BWA_VERSION=5961611c358

all: compile_assembler compile_index compile_assignment

clean:
	$(MAKE) -C prophyle_assembler clean
	$(MAKE) -C prophyle_index clean
	$(MAKE) -C prophyle_assignment clean

compile_assembler:
	$(MAKE) -C prophyle_assembler

compile_index:  prophyle_index/bwa/Makefile
	$(MAKE) -C prophyle_index

compile_assignment:
	$(MAKE) -C prophyle_assignment

submodules: prophyle_index/bwa/Makefile

prophyle_index/bwa/Makefile:
	cd prophyle_index/bwa && curl -L https://github.com/lh3/bwa/archive/$(BWA_VERSION).tar.gz | tar xz --strip-components 1
