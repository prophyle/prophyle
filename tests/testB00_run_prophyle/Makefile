.PHONY: all clean

include ../conf.mk

LIBRARIES = bacteria viruses plasmids hmp
SEQS = seqs
PKG_SEQS = pkg_seqs

$(LIBRARIES):
	@echo
	@echo "========================================="
	@echo "$@"
	@echo "========================================="
	@echo
	@echo "Download"
	@echo "--------"
	@echo
	$(PROP) download -F -d $(SEQS) -l $(SEQS)/log_down_$@.txt $@
	$(PROP_PKG) download -F -d $(PKG_SEQS) -l $(PKG_SEQS)/log_pkg_down_$@.txt $@
	@echo
	@echo "Index"
	@echo "-----"
	@echo
	$(PROP) index -F $(SEQS)/$@.nw idx_$@ 2>log_idx_$@.txt
	$(PROP_PKG) index -F $(PKG_SEQS)/$@.nw pkg_idx_$@ 2>log_idx_$@.txt
	@echo
	@echo "Classify"
	@echo "--------"
	@echo
	$(PROP) classify -P -A -f sam -m h1 idx_$@ $(FQ) $(FQ) >idx_$@/_test_class.sam
	$(PROP) classify -L -R -f kraken -m c1 idx_$@ $(FQ) >idx_$@/_test_class.kra
	$(PROP_PKG) classify -P -A -f sam -m h1 pkg_idx_$@ $(FQ) $(FQ) >pkg_idx_$@/_test_class.sam
	$(PROP_PKG) classify -L -R -f kraken -m c1 pkg_idx_$@ $(FQ) >pkg_idx_$@/_test_class.kra
	diff idx_$@/_test_class.sam pkg_idx_$@/_test_class.sam > diff_class_sam_$@.txt
	diff idx_$@/_test_class.kra pkg_idx_$@/_test_class.kra > diff_class_kra_$@.txt
	# ANALYZE

# Add LIBRARIES dependency to avoid parallel usage of the same seq directory
all_libs: $(LIBRARIES)
	@echo
	@echo "========================================="
	@echo "ALL"
	@echo "========================================="
	@echo
	@echo "Download"
	@echo "--------"
	@echo
	$(PROP) download -F -d $(SEQS) -l $(SEQS)/log_down_all.txt all
	$(PROP_PKG) download -F -d $(PKG_SEQS) -l $(PKG_SEQS)/log_pkg_down_all.txt all
	@echo
	@echo "Index"
	@echo "-----"
	@echo
	$(PROP) index -F -K -M -k 13 -s 0.05 $(SEQS)/*.nw idx_all 2>log_idx_all.txt
	$(PROP_PKG) index -F -K -M -k 13 -s 0.05 $(PKG_SEQS)/*.nw pkg_idx_all 2>log_pkg_idx_all.txt
	@echo
	@echo "Classify"
	@echo "--------"
	@echo
	$(PROP) classify idx_all $(FQ) >idx_all/_test_class.sam
	$(PROP_PKG) classify pkg_idx_all $(FQ) >pkg_idx_all/_test_class.sam
	diff idx_all/_test_class.sam pkg_idx_all/_test_class.sam > diff_class_all.txt
	# ANALYZE

git_pkg_diffs: $(LIBRARIES) all_libs
	@for f in ./diff_*.txt ; do test `wc -c < $$f` -eq 0 || (echo "file $$f is not empty" && exit 1) ; done

all: $(LIBRARIES) all_libs git_pkg_diffs

clean:
	rm -fr log_* idx_* pkg_idx_* $(SEQS) $(PKG_SEQS)
