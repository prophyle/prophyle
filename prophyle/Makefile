.PHONY: all clean compile_assembler compile_index

all: compile_assembler compile_index

clean:
	$(MAKE) -C prophyle_assembler clean
	$(MAKE) -C prophyle-index clean

compile_assembler:
	$(MAKE) -C prophyle_assembler

compile_index:
	$(MAKE) -C prophyle-index
