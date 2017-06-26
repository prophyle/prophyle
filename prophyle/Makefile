.PHONY: all clean compile_assembler compile_index

all: compile_assembler compile_index compile_assignment

clean:
	$(MAKE) -C prophyle_assembler clean
	$(MAKE) -C prophyle_index clean
	$(MAKE) -C prophyle_assignment clean

compile_assembler:
	$(MAKE) -C prophyle_assembler

compile_index:
	$(MAKE) -C prophyle_index

compile_assignment:
	$(MAKE) -C prophyle_assignment
