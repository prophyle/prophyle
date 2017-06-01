include params.mk


.PHONY: all clean

SHELL=/usr/bin/env bash
.SHELLFLAGS = -ecuf -o pipefail


PRG_ASM?=prophyle_assembler
PRG_DUST?=dustmasker

$(info )
$(info /------------------------------------------------------------------)

ifdef K
   $(info | K-mer length:           $(K))
else
   $(error | K-mer length is not specified)
endif

$(info | Assembler:              $(PRG_ASM))

$(info | DustMasker:             $(PRG_DUST))

ifdef MASKREP
   $(info | Masking repeats:        On)
   CMD_MASKING= | $(PRG_DUST) -infmt fasta -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g'
else
   $(info | Masking repeats:        Off)
   CMD_MASKING=
endif

ifdef NONPROP
   $(info | K-mer propagation:      Off)
else
   $(info | K-mer propagation:      On)
endif

ifdef NONDEL
   $(info | K-mer propagation mode: Non-deletative)
   REASM=1
else
   $(info | K-mer propagation mode: Deletative)
endif

ifdef REASM
   $(info | Re-assembling leaves:   On)
   CMD_REASM= | $(PRG_ASM) -S -i - -o -
else
   $(info | Re-assembling leaves:   Off)
   CMD_REASM=
endif
$(info \------------------------------------------------------------------)
$(info )

all: ./merge_root.reduced.fa.complete

clean:
	rm -f *.complete *.fa *.tsv

./merge_root.reduced.fa.complete: ./merge_root.full.fa.complete
	ln -s ./merge_root.full.fa ./merge_root.reduced.fa
	touch $@




./n1.full.fa.complete: /home/simone/git/prophyle/tests/testB02_check_results/fastas/n1.fa
	cat $^ $(CMD_MASKING) $(CMD_REASM) > ./n1.full.fa
	touch $@



./n2.full.fa.complete: /home/simone/git/prophyle/tests/testB02_check_results/fastas/n2.fa
	cat $^ $(CMD_MASKING) $(CMD_REASM) > ./n2.full.fa
	touch $@



./n3.full.fa.complete: /home/simone/git/prophyle/tests/testB02_check_results/fastas/n3.fa
	cat $^ $(CMD_MASKING) $(CMD_REASM) > ./n3.full.fa
	touch $@


ifdef NONDEL
   CMD_ASM_OUT_./n4.full.fa =
else
   CMD_ASM_OUT_./n4.full.fa = -o ./n1.reduced.fa -o ./n2.reduced.fa -o ./n3.reduced.fa
endif

ifdef NONPROP
   CMD_ASM_./n4.full.fa = touch ./n4.full.fa ./n1.reduced.fa ./n2.reduced.fa ./n3.reduced.fa
else
   CMD_ASM_./n4.full.fa = $(PRG_ASM) -S -k $(K) -i ./n1.full.fa -i ./n2.full.fa -i ./n3.full.fa $(CMD_ASM_OUT_./n4.full.fa) -x ./n4.full.fa -s ./n4.count.tsv
endif

./n4.full.fa.complete: ./n1.full.fa.complete ./n2.full.fa.complete ./n3.full.fa.complete
	@echo starting propagation for $@
	$(CMD_ASM_./n4.full.fa)
	touch $@


./n5.full.fa.complete: /home/simone/git/prophyle/tests/testB02_check_results/fastas/n5.fa
	cat $^ $(CMD_MASKING) $(CMD_REASM) > ./n5.full.fa
	touch $@



./n6.full.fa.complete: /home/simone/git/prophyle/tests/testB02_check_results/fastas/n6.fa
	cat $^ $(CMD_MASKING) $(CMD_REASM) > ./n6.full.fa
	touch $@


ifdef NONDEL
   CMD_ASM_OUT_./n7.full.fa =
else
   CMD_ASM_OUT_./n7.full.fa = -o ./n5.reduced.fa -o ./n6.reduced.fa
endif

ifdef NONPROP
   CMD_ASM_./n7.full.fa = touch ./n7.full.fa ./n5.reduced.fa ./n6.reduced.fa
else
   CMD_ASM_./n7.full.fa = $(PRG_ASM) -S -k $(K) -i ./n5.full.fa -i ./n6.full.fa $(CMD_ASM_OUT_./n7.full.fa) -x ./n7.full.fa -s ./n7.count.tsv
endif

./n7.full.fa.complete: ./n5.full.fa.complete ./n6.full.fa.complete
	@echo starting propagation for $@
	$(CMD_ASM_./n7.full.fa)
	touch $@

ifdef NONDEL
   CMD_ASM_OUT_./n8.full.fa =
else
   CMD_ASM_OUT_./n8.full.fa = -o ./n4.reduced.fa -o ./n7.reduced.fa
endif

ifdef NONPROP
   CMD_ASM_./n8.full.fa = touch ./n8.full.fa ./n4.reduced.fa ./n7.reduced.fa
else
   CMD_ASM_./n8.full.fa = $(PRG_ASM) -S -k $(K) -i ./n4.full.fa -i ./n7.full.fa $(CMD_ASM_OUT_./n8.full.fa) -x ./n8.full.fa -s ./n8.count.tsv
endif

./n8.full.fa.complete: ./n4.full.fa.complete ./n7.full.fa.complete
	@echo starting propagation for $@
	$(CMD_ASM_./n8.full.fa)
	touch $@

ifdef NONDEL
   CMD_ASM_OUT_./merge_root.full.fa =
else
   CMD_ASM_OUT_./merge_root.full.fa = -o ./n8.reduced.fa
endif

ifdef NONPROP
   CMD_ASM_./merge_root.full.fa = touch ./merge_root.full.fa ./n8.reduced.fa
else
   CMD_ASM_./merge_root.full.fa = $(PRG_ASM) -S -k $(K) -i ./n8.full.fa $(CMD_ASM_OUT_./merge_root.full.fa) -x ./merge_root.full.fa -s ./merge_root.count.tsv
endif

./merge_root.full.fa.complete: ./n8.full.fa.complete
	@echo starting propagation for $@
	$(CMD_ASM_./merge_root.full.fa)
	touch $@

