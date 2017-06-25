CXX      = g++
CXXFLAGS  = -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -g -O2

.PHONY: all clean 

all: prophyle_assignment

prophyle_assignment: prophyle_assignment.o
	$(CXX) $(CXXFLAGS) $(DFLAGS) prophyle_assignment.o read_processor.o tree_index.o knhx.o -o prophyle_assignment 

prophyle_assignment.o: read_processor.o tree_index.o knhx.o prophyle_assignment.cpp
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c read_processor.o knhx.o tree_index.o prophyle_assignment.cpp

read_processor.o: tree_index.o knhx.o read_processor.cpp read_processor.h
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c tree_index.o knhx.o read_processor.cpp read_processor.h

tree_index.o: knhx.o tree_index.cpp tree_index.h
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c knhx.o tree_index.cpp tree_index.h

knhx.o: knhx.c knhx.h
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c knhx.c knhx.h

clean:
	rm -f prophyle_assignment *.o
