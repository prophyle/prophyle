SHELL=/bin/bash -o pipefail

PROP_DIR=../../prophyle

ASM=$(PROP_DIR)/prophyle-assembler/prophyle-assembler
F2K=$(PROP_DIR)/_fa_to_kmers.py -m a
TEST_NEWICK=$(PROP_DIR)/test_newick_tree.py
