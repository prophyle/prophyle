SHELL=/bin/bash -o pipefail

PROP_DIR=../../prophyle

ASM=$(PROP_DIR)/prophyle-assembler/prophyle-assembler
IND=$(PROP_DIR)/prophyle-index/prophyle-index
F2K=$(PROP_DIR)/_fa_to_kmers.py -m a
1STEP=$(PROP_DIR)/1step_match.py
TEST_NEWICK=$(PROP_DIR)/test_newick_tree.py
