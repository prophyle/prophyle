#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "bntseq.h"
#include "bwa.h"
#include "bwt.h"
#include "utils.h"
#include "rle.h"
#include "rope.h"

#ifdef _DIVBWT
#include "divsufsort.h"
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is);
int bwa_pac2bwt(int argc, char *argv[]);
void bwt_bwtupdate_core(bwt_t *bwt);
int bwa_bwtupdate(int argc, char *argv[]);
int bwa_bwt2sa(int argc, char *argv[]);

#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

void bwt_revert_bwtupdate_core(bwt_t *bwt)
{
	bwtint_t i, k, c[4], n_occ;
	uint32_t *buf;

	bwt->bwt_size = (bwt->seq_len + 15) >> 4;
	// n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	// bwt->bwt_size += n_occ * sizeof(bwtint_t); // the new size
	buf = (uint32_t*)calloc(bwt->bwt_size, 4); // will be the new bwt
	for (i = k = 0; i < bwt->seq_len; ++i) {
		if (i % OCC_INTERVAL == 0) {
			k += sizeof(bwtint_t); // in fact: sizeof(bwtint_t)=4*(sizeof(bwtint_t)/4)
		}
		if (i % 16 == 0) {
			buf[i / 16] = bwt->bwt[k]; // 16 == sizeof(uint32_t)/2
			k++;
		}
	}
	// the last element
	// update bwt
	free(bwt->bwt); bwt->bwt = buf;
}

int bwa_build_compressed_index(const char *fa, const char *prefix, int algo_type, int block_size)
{
	char *str, *str2;
	clock_t t;
	int64_t l_pac;

	str  = (char*)calloc(strlen(prefix) + 10, 1);
	str2 = (char*)calloc(strlen(prefix) + 10, 1);

	// { // nucleotide indexing
	// 	gzFile fp = xzopen(fa, "r");
	// 	t = clock();
	// 	if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Pack FASTA... ");
	// 	l_pac = bns_fasta2bntseq(fp, prefix, 0);
	// 	if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	// 	err_gzclose(fp);
	// }
	//if (algo_type == 0) algo_type = l_pac > 50000000? 2 : 3; // set the algorithm for generating BWT
	// {
	// 	strcpy(str, prefix); strcat(str, ".pac");
	// 	strcpy(str2, prefix); strcat(str2, ".bwt");
	// 	t = clock();
	// 	if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Construct BWT for the packed sequence...\n");
	// 	if (algo_type == 2) bwt_bwtgen2(str, str2, block_size);
	// 	else if (algo_type == 1 || algo_type == 3) {
	// 		bwt_t *bwt;
	// 		bwt = bwt_pac2bwt(str, algo_type == 3);
	// 		fprintf(stderr, "bwt->seq_len = %d\n", bwt->seq_len);
	// 		fprintf(stderr, "bwt->bwt_size = %d\n", bwt->bwt_size);
	// 		for (int i = 0; i < bwt->bwt_size && i < 10; ++i) {
	// 			fprintf(stderr, "%u\n", bwt->bwt[i]);
	// 		}
	// 		bwt_dump_bwt(str2, bwt);
	// 		bwt_destroy(bwt);
	// 	}
	// 	if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	// }
	{
    bwt_t *bwt;
    strcpy(str, prefix); strcat(str, ".bwt");
    t = clock();
    if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Revert update BWT... ");
    bwt = bwt_restore_bwt(str);
    bwt_revert_bwtupdate_core(bwt);
    bwt_dump_bwt(str, bwt);
    bwt_destroy(bwt);
    if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  }
	free(str); free(str2);
	return 0;
}

int bwa_compress(int argc, char *argv[]) // the "compress" command
{
	int c, algo_type = BWTALGO_AUTO, block_size = 10000000;
	char *prefix = 0, *str;
	while ((c = getopt(argc, argv, "a:p:b:")) >= 0) {
		switch (c) {
		case 'a': // if -a is not set, algo_type will be determined later
			if (strcmp(optarg, "rb2") == 0) algo_type = BWTALGO_RB2;
			else if (strcmp(optarg, "bwtsw") == 0) algo_type = BWTALGO_BWTSW;
			else if (strcmp(optarg, "is") == 0) algo_type = BWTALGO_IS;
			else err_fatal(__func__, "unknown algorithm: '%s'.", optarg);
			break;
		case 'p': prefix = strdup(optarg); break;
		case 'b':
			block_size = strtol(optarg, &str, 10);
			if (*str == 'G' || *str == 'g') block_size *= 1024 * 1024 * 1024;
			else if (*str == 'M' || *str == 'm') block_size *= 1024 * 1024;
			else if (*str == 'K' || *str == 'k') block_size *= 1024;
			break;
		default: return 1;
		}
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa compress [options] <in.fasta>\n\n");
		fprintf(stderr, "Options: -a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]\n");
		fprintf(stderr, "         -p STR    prefix of the index [same as fasta name]\n");
		fprintf(stderr, "         -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [%d]\n", block_size);
		fprintf(stderr, "         -6        index files named as <in.fasta>.64.* instead of <in.fasta>.* \n");
		fprintf(stderr, "\n");
		fprintf(stderr,	"Warning: `-a bwtsw' does not work for short genomes, while `-a is' and\n");
		fprintf(stderr, "         `-a div' do not work not for long genomes.\n\n");
		return 1;
	}
	if (prefix == 0) {
		prefix = malloc(strlen(argv[optind]) + 4);
		strcpy(prefix, argv[optind]);
	}
	bwa_build_compressed_index(argv[optind], prefix, algo_type, block_size);
	free(prefix);
	return 0;
}

#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

int bwa_decompress_index(const char *fa, const char *prefix) {
  char *str, *str3;
  clock_t t;
  str  = (char*)calloc(strlen(prefix) + 10, 1);
  str3  = (char*)calloc(strlen(prefix) + 10, 1);

  {
    bwt_t *bwt;
    strcpy(str, prefix); strcat(str, ".bwt");
    t = clock();
    if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Update BWT... ");
    bwt = bwt_restore_bwt(str);
    bwt_bwtupdate_core(bwt);
    bwt_dump_bwt(str, bwt);
    bwt_destroy(bwt);
    if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  }
  {
    gzFile fp = xzopen(fa, "r");
    t = clock();
    if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Pack forward-only FASTA... ");
    bns_fasta2bntseq(fp, prefix, 1);
    if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    err_gzclose(fp);
  }
  {
    bwt_t *bwt;
    strcpy(str, prefix); strcat(str, ".bwt");
    strcpy(str3, prefix); strcat(str3, ".sa");
    t = clock();
    if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Construct SA from BWT and Occ... ");
    bwt = bwt_restore_bwt(str);
    bwt_cal_sa(bwt, 32);
    bwt_dump_sa(str3, bwt);
    bwt_destroy(bwt);
    if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  }
  free(str); free(str3);
  return 0;
}

int bwa_decompress(int argc, char *argv[]) // the "decompress" command
{
	int c;
	char *prefix = 0;
	while ((c = getopt(argc, argv, "p:")) >= 0) {
		switch (c) {
		case 'p': prefix = strdup(optarg); break;
		default: return 1;
		}
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa decompress [options] <in.fasta>\n\n");
		fprintf(stderr, "Options: -p STR    prefix of the index [same as fasta name]\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (prefix == 0) {
		prefix = malloc(strlen(argv[optind]) + 4);
		strcpy(prefix, argv[optind]);
	}
	bwa_decompress_index(argv[optind], prefix);
	free(prefix);
	return 0;
}
