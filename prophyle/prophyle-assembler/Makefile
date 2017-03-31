CXX      ?= g++
CXXFLAGS  = -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -g -O2
LIBS      = -lz

.PHONY: all clean 

all: prophyle-assembler

prophyle-assembler: prophyle-assembler.o
	$(CXX) $(CXXFLAGS) $(DFLAGS) $^ -o $@ -L. $(LIBS)

prophyle-assembler.o: prophyle-assembler.cpp kseq.h
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $<

clean:
	rm -f prophyle-assembler *.o
