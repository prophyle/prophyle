CXX      ?= g++
CXXFLAGS  = -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -g -O2
LIBS      = -lz

.PHONY: all clean assembler

all: assembler

assembler: main.o
	$(CXX) $(CXXFLAGS) $(DFLAGS) $^ -o $@ -L. $(LIBS)

main.o: main.cpp kseq.h
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $<

clean:
	rm -f assembler *.o