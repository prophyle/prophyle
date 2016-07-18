.PHONY: all test utils clean

CXX=	g++
RM=		rm -f

SRCS=	krakendb.cpp quickfile.cpp
OBJS=	$(subst .cpp,.o,$(SRCS))
TEST_DIR = test
TEST_DB= test_db.kdb
TEST_SORTED= $(TEST_DIR)/test_sorted.kdb
TEST_SORTED_IDX = $(TEST_DIR)/test_sorted.idx
TEST_TAXA_MAP = $(TEST_DIR)/test_taxa_map.dmp

all: get_kmers_by_taxid taxid_sorter

get_kmers_by_taxid: $(OBJS)
	$(CXX) -o $@ get_kmers_by_taxid.cpp $(OBJS)

taxid_sorter: $(OBJS)
	$(CXX) -o $@ taxid_sorter.cpp $(OBJS)

krakendb.o: krakendb.hpp quickfile.hpp kraken_headers.hpp

quickfile.o: quickfile.hpp kraken_headers.hpp

utils: taxa_list

taxa_list: taxa_list.cpp kraken_headers.hpp
	$(CXX) -o $@ $<

test: all utils ./test.sh
	mkdir -p $(TEST_DIR) && \
	./taxid_sorter -d $(TEST_DB) -o $(TEST_SORTED) \
		-i $(TEST_SORTED_IDX) -m $(TEST_TAXA_MAP) && \
	./taxa_list -m $(TEST_TAXA_MAP) -o $(TEST_TAXA_LIST) && \
	./test.sh $(TEST_DB) $(TEST_SORTED) $(TEST_SORTED_IDX) \
		$(TEST_TAXA_MAP) $(TEST_TAXA_LIST)

clean:
	$(RM) $(OBJS) *~
