.PHONY: all clean

include ../conf.mk

K=2
tree=tree.nw
index=index

all:
	$(PROP) index -A -g fastas -k 2 $(tree) $(index)

clean:
	rm -rf $(index) $(res) $(diff_res)
