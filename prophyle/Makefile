.PHONY: all clean compile_assembler compile_index

all: compile_assembler compile_index fake_extension.c

clean:
	$(MAKE) -C prophyle-assembler clean
	$(MAKE) -C prophyle-index clean

compile_assembler:
	$(MAKE) -C prophyle-assembler

compile_index:
	$(MAKE) -C prophyle-index


fake_extension.c:
	echo "int main () {}" > $@