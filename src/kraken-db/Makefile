.PHONY: all clean

CXX=	g++
RM=		rm -f

SRCS=	krakendb.cpp quickfile.cpp
OBJS=	$(subst .cpp,.o,$(SRCS))

all: get_kmers_by_taxid taxid_sorter

get_kmers_by_taxid: $(OBJS)
	$(CXX) -o $@ get_kmers_by_taxid.cpp $(OBJS)

taxid_sorter: $(OBJS)
	$(CXX) -o $@ taxid_sorter.cpp $(OBJS)

krakendb.o: krakendb.hpp quickfile.hpp kraken_headers.hpp

quickfile.o: quickfile.hpp kraken_headers.hpp

clean:
	$(RM) $(OBJS) *~
