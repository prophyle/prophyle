#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "prophyle_query.h"
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"
#include "bwa.h"
#include "bwase.h"
#include "kstring.h"
#include "klcp.h"
#include "bwa_utils.h"
#include "contig_node_translator.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.1"
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

int exk_match(int argc, char *argv[])
{
	int c, opte = -1;
	exk_opt_t *opt;
	char *prefix;

	opt = exk_init_opt();
	while ((c = getopt(argc, argv, "l:psuvk:bt:")) >= 0) {
		switch (c) {
		case 'v': { opt->output_old = 1; opt->output = 0; } break;
		case 'u': opt->use_klcp = 1; break;
		case 'k': opt->kmer_length = atoi(optarg); break;
		case 's': opt->skip_after_fail = 1; break;
		case 'p': opt->skip_positions_on_border = 0; break;
		case 'l': { opt->need_log = 1; opt->log_file_name = optarg; break; }
		case 'b': opt->output_read_qual = 1; break;
		case 't': opt->n_threads = atoi(optarg); break;
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (opt->output_old && opt->n_threads > 1) {
		fprintf(stderr, "-v option can be used only with one thread (-t 1)\n");
		return 1;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   prophyle-index query [options] <prefix> <in.fq>\n\n");
		fprintf(stderr, "Options: -k INT    length of k-mer\n");
		fprintf(stderr, "         -u        use klcp for matching\n");
		fprintf(stderr, "         -v        output set of chromosomes for every k-mer\n");
		fprintf(stderr, "         -p        do not check whether k-mer is on border of two contigs, and show such k-mers in output\n");
		fprintf(stderr, "         -l char*  log file name to output statistics\n");
		fprintf(stderr, "         -t INT    number of threads [%d]\n", opt->n_threads);
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

int exk_index(int argc, char *argv[])
{
	int c, opte = -1;	exk_opt_t *opt;
	char *prefix;
	opt = exk_init_opt();
	int sa_intv = 32;
	while ((c = getopt(argc, argv, "si:k:")) >= 0) {
		switch (c) {
		case 'k': opt->kmer_length = atoi(optarg); break;
		case 'i': sa_intv = atoi(optarg); break;
		case 's': opt->construct_sa_parallel = 1; break;
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   prophyle-index build <prefix>\n\n");
		// fprintf(stderr, "Options: -t INT    number of threads [%d]\n", opt->n_threads);
		// fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
    // fprintf(stderr, "         -f FILE   file to write output to instead of stdout\n");

		fprintf(stderr, "Options:  -k INT    length of k-mer\n");
		fprintf(stderr, "          -s        construct klcp and sa in parallel\n");
		fprintf(stderr, "          -i        sampling distance for SA\n");
    fprintf(stderr, "\n");
		return 1;
	}
	if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
		fprintf(stderr, "[%s] fail to locate the index %s\n", __func__, argv[optind]);
		return 1;
	}
	exk_index_core(prefix, opt, sa_intv);
	free(prefix);
	return 0;
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: prophyle-index (alignment of k-mers)\n");
	fprintf(stderr, "Usage:   prophyle-index command [options]\n\n");
	fprintf(stderr, "Command: build         construct index\n");
	fprintf(stderr, "Command: query         query reads against index\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	extern char *bwa_pg;
	int i, ret = 0;
	kstring_t pg = {0,0,0};
	ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	bwa_pg = pg.s;
	if (argc < 2) return usage();
	if (strcmp(argv[1], "build") == 0) ret = exk_index(argc - 1, argv + 1);
	else if (strcmp(argv[1], "query") == 0) ret = exk_match(argc-1, argv+1);
	else return usage();

	err_fflush(stdout);
	err_fclose(stdout);

	free(bwa_pg);
	return ret;
}
