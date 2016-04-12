#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtexk.h"
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"
#include "bwa.h"
#include "bwase.h"
#include "kstring.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.1"
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

exk_opt_t *exk_init_opt()
{
	exk_opt_t *o;
	o = (exk_opt_t*)calloc(1, sizeof(exk_opt_t));
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->n_threads = 1;
	o->trim_qual = 0;
	o->kmer_length = 14;
	return o;
}

int bwt_cal_sa_coord(const bwt_t *bwt, int len, const ubyte_t *str, int* k, int* l, int start_pos)
{
	bwtint_t ok, ol;
	int i, bid;
	bid = 0;
	*k = 0; *l = bwt->seq_len;

	//fprintf(stderr, "start k = %d, l = %d\n", *k, *l);
	for (i = start_pos; i < start_pos + len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(bwt, *k - 1, *l, c, &ok, &ol);
			*k = bwt->L2[c] + ok + 1;
			*l = bwt->L2[c] + ol;
			//fprintf(stderr, "after i = %d character, cur k = %d, l = %d\n", i, *k, *l);
		}
		if (*k > *l || c > 3) { // then restart
			return i - start_pos;
		}
	}
	return len;
}

void bwa_cal_sa(int tid, bwaidx_t* idx, int n_seqs, bwa_seq_t *seqs, const exk_opt_t *opt)
{
	bwase_initialize();
	int i, j;

	bwt_t* bwt = idx->bwt;

	fprintf(stdout, "\n");
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
		p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;

		// NEED TO UNDERSTAND
		// core function
		// for (j = 0; j < p->len; ++j) // we need to complement
		// 	p->seq[j] = p->seq[j] > 3? 4 : 3 - p->seq[j];

		//fprintf(stderr, "seq = ");
		fprintf(stdout, "#");
		for(j = 0; j < p->len; ++j) {
			fprintf(stdout, "%c", "ACGTN"[p->seq[j]]);
		}
		fprintf(stdout, "\n");
		int k, l;
		for(int start_pos = p->len - opt->kmer_length; start_pos >= 0; --start_pos) {
			k = 0;
			l = 0;
			bwt_cal_sa_coord(bwt, opt->kmer_length, p->seq, &k, &l, start_pos);
			//fprintf(stderr, "start_pos = %d\n", start_pos);
			//fprintf(stderr, "found k = %d, l = %d\n", k, l);
			int* seen_rids = malloc((l - k + 1) * sizeof(int));
			int rids_cnt = 0;
			for(int t = k; t <= l; ++t) {
				int strand;
				int pos = bwa_sa2pos(idx->bns, idx->bwt, t, p->len, &strand);//bwt_sa(bwt, t);
				int rid = bns_pos2rid(idx->bns, pos);
				int seen = 0;
				for(int r = 0; r < rids_cnt; ++r) {
					if (seen_rids[r] == rid) {
						seen = 1;
						break;
					}
				}
				if (!seen && rid != -1) {
					seen_rids[rids_cnt] = rid;
					++rids_cnt;
					//fprintf(stderr, "t = %d, pos = %d, rid = %d\n", t, pos, rid);
				}
			}
			fprintf(stdout, "%d ", rids_cnt);
			for(int r = 0; r < rids_cnt; ++r) {
				fprintf(stdout, "%d ", seen_rids[r]);
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "#\n");


		free(p->name); free(p->seq); free(p->rseq); free(p->qual);
		p->name = 0; p->seq = p->rseq = p->qual = 0;
	}
}

bwa_seqio_t *bwa_open_reads_new(int mode, const char *fn_fa)
{
	bwa_seqio_t *ks;
	if (mode & BWA_MODE_BAM) { // open BAM
		int which = 0;
		if (mode & BWA_MODE_BAM_SE) which |= 4;
		if (mode & BWA_MODE_BAM_READ1) which |= 1;
		if (mode & BWA_MODE_BAM_READ2) which |= 2;
		if (which == 0) which = 7; // then read all reads
		ks = bwa_bam_open(fn_fa, which);
	} else ks = bwa_seq_open(fn_fa);
	return ks;
}

void bwa_exk_core(const char *prefix, const char *fn_fa, const exk_opt_t *opt) {
	int n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	bwt_t *bwt;
	bwaidx_t* idx;

	// initialization
	{ // load BWT
		//fprintf(stderr, "%s\n", prefix);
		if ((idx = bwa_idx_load(prefix, BWA_IDX_ALL)) == 0) {
			fprintf(stderr, "Couldn't load idx from %s\n", prefix);
			return;
		}
		bwt = idx->bwt;
	}

	ks = bwa_open_reads_new(opt->mode, fn_fa);
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
		tot_seqs += n_seqs;
		t = clock();

		//fprintf(stderr, "[bwa_aln_core] calculate SA coordinate...\n ");

		bwa_cal_sa(0, idx, n_seqs, seqs, opt);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		bwa_free_read_seq(n_seqs, seqs);
		//fprintf(stderr, "[bwa_aln_core] %d sequences have been processed.\n", tot_seqs);
	}

	// destroy
	bwt_destroy(bwt);
	bwa_seq_close(ks);
}


int bwa_exk(int argc, char *argv[])
{
	int c, opte = -1;
	exk_opt_t *opt;
	char *prefix;

	opt = exk_init_opt();
	while ((c = getopt(argc, argv, "a:n:o:e:i:d:l:k:LR:m:t:NM:O:E:q:f:b012IYB:")) >= 0) {
		switch (c) {
		case 'a': opt->kmer_length = atoi(optarg); break;
		case 'e': opte = atoi(optarg); break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
		case 'q': opt->trim_qual = atoi(optarg); break;
		case 'N': opt->mode |= BWA_MODE_NONSTOP; break;
		case 'f': xreopen(optarg, "wb", stdout); break;
		case 'b': opt->mode |= BWA_MODE_BAM; break;
		case '0': opt->mode |= BWA_MODE_BAM_SE; break;
		case '1': opt->mode |= BWA_MODE_BAM_READ1; break;
		case '2': opt->mode |= BWA_MODE_BAM_READ2; break;
		case 'I': opt->mode |= BWA_MODE_IL13; break;
		case 'Y': opt->mode |= BWA_MODE_CFY; break;
		case 'B': opt->mode |= atoi(optarg) << 24; break;
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   exk exk [options] <prefix> <in.fq>\n\n");
		fprintf(stderr, "Options: -t INT    number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
    fprintf(stderr, "         -f FILE   file to write output to instead of stdout\n");
		fprintf(stderr, "         -B INT    length of barcode\n");
		fprintf(stderr, "         -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
		fprintf(stderr, "         -b        the input read file is in the BAM format\n");
		fprintf(stderr, "         -0        use single-end reads only (effective with -b)\n");
		fprintf(stderr, "         -1        use the 1st read in a pair (effective with -b)\n");
		fprintf(stderr, "         -2        use the 2nd read in a pair (effective with -b)\n");
		fprintf(stderr, "         -Y        filter Casava-filtered sequences\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
		fprintf(stderr, "[%s] fail to locate the index %s\n", __func__, argv[optind]);
		free(opt);
		return 1;
	}
	bwa_exk_core(prefix, argv[optind+1], opt);
	free(opt); free(prefix);
	return 0;
}


static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: bwa (alignment via Burrows-Wheeler transformation)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Heng Li <lh3@sanger.ac.uk>\n\n");
	fprintf(stderr, "Usage:   bwa <command> [options]\n\n");
	fprintf(stderr, "Command: index         index sequences in the FASTA format\n");
	fprintf(stderr, "         mem           BWA-MEM algorithm\n");
	fprintf(stderr, "         fastmap       identify super-maximal exact matches\n");
	fprintf(stderr, "         pemerge       merge overlapping paired ends (EXPERIMENTAL)\n");
	fprintf(stderr, "         aln           gapped/ungapped alignment\n");
	fprintf(stderr, "         samse         generate alignment (single ended)\n");
	fprintf(stderr, "         sampe         generate alignment (paired ended)\n");
	fprintf(stderr, "         bwasw         BWA-SW for long queries\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "         shm           manage indices in shared memory\n");
	fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
	fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
	fprintf(stderr, "         pac2bwtgen    alternative algorithm for generating BWT\n");
	fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
	fprintf(stderr, "         bwt2sa        generate SA from BWT and Occ\n");
	fprintf(stderr, "\n");
	fprintf(stderr,
"Note: To use BWA, you need to first index the genome with `bwa index'.\n"
"      There are three alignment algorithms in BWA: `mem', `bwasw', and\n"
"      `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'\n"
"      first. Please `man ./bwa.1' for the manual.\n\n");
	return 1;
}

int main(int argc, char *argv[])
{
	extern char *bwa_pg;
	int i, ret;
	double t_real;
	kstring_t pg = {0,0,0};
	t_real = realtime();
	ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	bwa_pg = pg.s;
	if (argc < 2) return usage();
	else ret = bwa_exk(argc, argv);

	err_fflush(stdout);
	err_fclose(stdout);
	if (ret == 0) {
		//fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		//fprintf(stderr, "[%s] CMD:", __func__);
		//for (i = 0; i < argc; ++i)
		//	fprintf(stderr, " %s", argv[i]);
		//fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	}
	free(bwa_pg);
	return ret;
}
