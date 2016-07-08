CXX=	g++
RM=		rm -f

SRCS=	fast_converter.cpp krakendb.cpp quickfile.cpp
OBJS=	$(subst .cc,.o,$(SRCS))

all: converter

converter: $(OBJS)
	$(CXX) -o $@ $(OBJS)

converter.o: fast_converter.cpp krakendb.hpp quickfile.hpp kraken_headers.hpp

krakendb.o: krakendb.hpp quickfile.hpp kraken_headers.hpp

quickfile.o: quickfile.cpp kraken_headers.hpp

clean:
	$(RM) $(OBJS) *~
