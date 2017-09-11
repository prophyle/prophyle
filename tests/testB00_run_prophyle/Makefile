LIBRARIES = bacteria viruses plasmids hmp
GIT_LIBS = $(addprefix .complete.git.,$(LIBRARIES))
PKG_LIBS = $(addprefix .complete.pkg.,$(LIBRARIES))
GIT_SEQS = seqs
PKG_SEQS = pkg_seqs
STATS = w u wl ul

.PHONY: all git pkg clean diffs

include ../conf.mk

$(GIT_LIBS):
	$(eval CUR_LIB := $(subst .complete.git.,,$@))
	@echo
	@echo "========================================="
	@echo "$(CUR_LIB)-git"
	@echo "========================================="
	@echo
	@echo "Download"
	@echo "--------"
	@echo
	$(PROP) download -F -d $(GIT_SEQS) -l $(GIT_SEQS)/log_down_$(CUR_LIB).txt $(CUR_LIB)
	@echo
	@echo "Index"
	@echo "-----"
	@echo
	$(PROP) index -F $(GIT_SEQS)/$(CUR_LIB).nw idx_$(CUR_LIB) 2>log_idx_$(CUR_LIB).txt
	@echo
	@echo "Classify"
	@echo "--------"
	@echo
	$(PROP) classify -P -A -f sam -m h1 idx_$(CUR_LIB) $(FQ) $(FQ) >idx_$(CUR_LIB)/_test_class.sam
	$(PROP) classify -L -K -f kraken -m c1 idx_$(CUR_LIB) $(FQ) >idx_$(CUR_LIB)/_test_class.kra
	@echo
	@echo "Analyze"
	@echo "-------"
	@echo
	@for s in $(STATS) ; do (cd idx_$(CUR_LIB) && ../$(PROP) analyze -s $$s . _test_analyze_$$s _test_class.sam) ; done
	touch $@

$(PKG_LIBS):
	$(eval CUR_LIB := $(subst .complete.pkg.,,$@))
	@echo
	@echo "========================================="
	@echo "$(CUR_LIB)-pkg"
	@echo "========================================="
	@echo
	@echo "Download"
	@echo "--------"
	@echo
	$(PROP_PKG) download -F -d $(PKG_SEQS) -l $(PKG_SEQS)/log_pkg_down_$(CUR_LIB).txt $(CUR_LIB)
	@echo
	@echo "Index"
	@echo "-----"
	@echo
	$(PROP_PKG) index -F $(PKG_SEQS)/$(CUR_LIB).nw pkg_idx_$(CUR_LIB) 2>log_idx_$(CUR_LIB).txt
	@echo
	@echo "Classify"
	@echo "--------"
	@echo
	$(PROP_PKG) classify -P -A -f sam -m h1 pkg_idx_$(CUR_LIB) $(FQ) $(FQ) >pkg_idx_$(CUR_LIB)/_test_class.sam
	$(PROP_PKG) classify -L -K -f kraken -m c1 pkg_idx_$(CUR_LIB) $(FQ) >pkg_idx_$(CUR_LIB)/_test_class.kra
	@echo
	@echo "Analyze"
	@echo "-------"
	@echo
	@for s in $(STATS) ; do (cd pkg_idx_$(CUR_LIB) && $(PROP_PKG) analyze -s $$s . _test_analyze_$$s _test_class.sam) ; done
	touch $@

# Add GIT_LIBS dependency to avoid parallel usage of the same seq directory
.complete.git.all: $(GIT_LIBS)
	@echo
	@echo "========================================="
	@echo "ALL-GIT"
	@echo "========================================="
	@echo
	@echo "Download"
	@echo "--------"
	@echo
	$(PROP) download -d $(GIT_SEQS) -l $(GIT_SEQS)/log_down_all.txt all
	@echo
	@echo "Index"
	@echo "-----"
	@echo
	$(PROP) index -F -K -M -k 13 -s 0.05 $(GIT_SEQS)/*.nw idx_all 2>log_idx_all.txt
	@echo
	@echo "Classify"
	@echo "--------"
	@echo
	$(PROP) classify idx_all $(FQ) >idx_all/_test_class.sam
	@echo
	@echo "Analyze"
	@echo "-------"
	@echo
	@for s in $(STATS) ; do (cd idx_all && ../$(PROP) analyze -s $$s . _test_analyze_$$s _test_class.sam) ; done
	touch $@

# Add PKG_LIBS dependency to avoid parallel usage of the same seq directory
.complete.pkg.all: $(PKG_LIBS)
	@echo
	@echo "========================================="
	@echo "ALL-PKG"
	@echo "========================================="
	@echo
	@echo "Download"
	@echo "--------"
	@echo
	$(PROP_PKG) download -d $(PKG_SEQS) -l $(PKG_SEQS)/log_pkg_down_all.txt all
	@echo
	@echo "Index"
	@echo "-----"
	@echo
	$(PROP_PKG) index -F -K -M -k 13 -s 0.05 $(PKG_SEQS)/*.nw pkg_idx_all 2>log_pkg_idx_all.txt
	@echo
	@echo "Classify"
	@echo "--------"
	@echo
	$(PROP_PKG) classify pkg_idx_all $(FQ) >pkg_idx_all/_test_class.sam
	@echo
	@echo "Analyze"
	@echo "-------"
	@echo
	@for s in $(STATS) ; do (cd pkg_idx_all && $(PROP_PKG) analyze -s $$s . _test_analyze_$$s _test_class.sam) ; done
	touch $@

diffs: $(GIT_LIBS) $(PKG_LIBS) .complete.git.all .complete.pkg.all
	@for l in $(LIBRARIES) all ; do diff idx_$$l/_test_class.sam pkg_idx_$$l/_test_class.sam > diff_class_sam_$$l.txt ; done
	@for l in $(LIBRARIES) ; do diff idx_$$l/_test_class.kra pkg_idx_$$l/_test_class.kra > diff_class_kra_$$l.txt ; done
	@for f in ./diff_*.txt ; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

git: $(GIT_LIBS)

pkg: $(PKG_LIBS)

all: git pkg .complete.git.all .complete.pkg.all diffs

clean:
	rm -fr .complete.* log_* idx_* pkg_idx_* $(GIT_SEQS) $(PKG_SEQS)
