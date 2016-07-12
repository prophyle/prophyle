CXX=	g++
RM=		rm -f

SRCS=	krakendb_taxid_sorter.cpp krakendb.cpp quickfile.cpp
OBJS=	$(subst .cpp,.o,$(SRCS))

all: krakendb_taxid_sorter

krakendb_taxid_sorter: $(OBJS)
	$(CXX) -o $@ $(OBJS)

krakendb_taxid_sorter.o: krakendb_taxid_sorter.cpp krakendb.hpp quickfile.hpp kraken_headers.hpp

krakendb.o: krakendb.hpp quickfile.hpp kraken_headers.hpp

quickfile.o: quickfile.cpp kraken_headers.hpp

clean:
	$(RM) $(OBJS) *~
