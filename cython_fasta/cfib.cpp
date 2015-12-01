#include "cfib.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)
	
	int fib(char *fasta) {
		gzFile fp;
		kseq_t *seq;
		int l;
		fp = gzopen(fasta, "r");
		seq = kseq_init(fp);
		while ((l = kseq_read(seq)) >= 0) {
			//printf("name: %s\n", seq->name.s);
			if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
			//printf("seq: %s\n", seq->seq.s);
			if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
		}


		kseq_destroy(seq);
		gzclose(fp);

		return 5;
	}
