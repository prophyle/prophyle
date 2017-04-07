CXX      ?= g++
CXXFLAGS  = -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -g -O2
LIBS      = -lz

.PHONY: all clean 

all: prophyle_assembler

prophyle_assembler: prophyle_assembler.o
	$(CXX) $(CXXFLAGS) $(DFLAGS) $^ -o $@ -L. $(LIBS)

prophyle_assembler.o: prophyle_assembler.cpp kseq.h
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $<

clean:
	rm -f prophyle_assembler *.o
