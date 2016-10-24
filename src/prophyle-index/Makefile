CC=			gcc
#CC=			clang --analyze
CFLAGS=		-g -Wall -Wno-unused-function -O2
WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
AR=			ar
DFLAGS=		-DHAVE_PTHREAD $(WRAP_MALLOC)

PROG=prophyle-index

AOBJS2=	\
			bwa/bwashm.o \
			bwa/bwase.o \
			bwa/bwaseqio.o \
			bwa/bwtgap.o \
			bwa/bwtaln.o \
			bwa/bamlite.o \
			bwa/bwape.o \
			bwa/kopen.o \
			bwa/pemerge.o \
			bwa/maxk.o \
			bwa/bwtsw2_core.o \
			bwa/bwtsw2_main.o \
			bwa/bwtsw2_aux.o \
			bwa/bwt_lite.o \
			bwa/bwtsw2_chain.o \
			bwa/fastmap.o \
			bwa/bwtsw2_pair.o \


INCLUDES=	-Ibwa
LIBS=		-lm -lz -lpthread
SUBDIRS=	.

ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

clean:
	rm -f gmon.out *.o a.out $(PROG) *~ *.a
	$(MAKE) -C bwa clean

$(PROG): bwa/libbwa.a $(AOBJS2) main.o prophyle_query.o klcp.o bitarray.o bwa_utils.o contig_node_translator.o
	$(CC) $(INCLUDES) $(CFLAGS) $(DFLAGS) $(AOBJS2) main.o prophyle_query.o klcp.o bitarray.o bwa_utils.o contig_node_translator.o -o $@ -Lbwa -lbwa $(LIBS)

#bwa/libbwa.a $(AOBJS2) bwtexk.o:
bwa/libbwa.a:
	 $(MAKE) -C bwa
