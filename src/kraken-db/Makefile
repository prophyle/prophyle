.PHONY: all clean

CXX=	g++
RM=		rm -f

SRCS=	krakendb.cpp quickfile.cpp
OBJS=	$(subst .cpp,.o,$(SRCS))

all: get_kmers_by_taxid taxid_sorter taxa_list

get_kmers_by_taxid: $(OBJS)
	$(CXX) -o $@ $@.cpp $(OBJS)

taxid_sorter: $(OBJS)
	$(CXX) -o $@ $@.cpp $(OBJS)

taxa_list:
	$(CXX) -o $@ $@.cpp

$(OBJS): %.o : %.cpp
	$(CXX) -c -o $@ $<

clean:
	$(RM) $(OBJS) *~
