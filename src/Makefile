.PHONY: all clean compile_assembler compile_index #compile_ext

all: compile_assembler compile_index #compile_ext

clean:
	$(MAKE) -C prophyle-assembler clean
	$(MAKE) -C prophyle-index clean
	#$(MAKE) -C kraken-db clean

compile_assembler:
	$(MAKE) -C prophyle-assembler

compile_index:
	$(MAKE) -C prophyle-index

#compile_ext:
#	$(MAKE) -C ../ext
