.PHONY: all clean

CXX=	g++
RM=		rm -rf

SRCS=	krakendb.cpp quickfile.cpp
OBJS=	$(subst .cpp,.o,$(SRCS))
TEST=	test
TMP=	tmp
UTILS=	$(TEST)/build_test_db $(TEST)/fast_converter

all: get_kmers_by_taxid taxid_sorter

get_kmers_by_taxid: $(OBJS)
	$(CXX) -o $@ $@.cpp $(OBJS)

taxid_sorter: $(OBJS)
	$(CXX) -o $@ $@.cpp $(OBJS)

$(UTILS): % : %.cpp $(OBJS)
	$(CXX) -o $@ $< $(OBJS)

$(OBJS): %.o : %.cpp
	$(CXX) -c -o $@ $<

test: all $(UTILS)
	test/test.sh $(TMP)

clean:
	$(RM) $(OBJS) *~
	$(RM) $(TEST)/$(TMP)
