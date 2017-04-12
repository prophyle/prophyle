/*
	prophyle_index main file, implements build and query commands.
	Author: Kamil Salikhov <salikhov.kamil@gmail.com>
	Licence: MIT
	Examples:
		build prophyle index for k = 20, suffix array and klcp simultaneously:
			prophyle_index build -k 20 -s index.fa
		query reads for k=20 using rolling window search with 10 threads, writing output in results.txt:
			prophyle_index query -u -k 20 -t 10 index.fa reads.fq > results.txt
*/

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "prophyle_query.h"
#include "bwa.h"
#include "prophyle_index_build.h"
#include "bwa_utils.h"

int prophyle_index_query(int argc, char *argv[])
{
	int c;
	prophyle_index_opt_t *opt;
	char *prefix;

	opt = prophyle_index_init_opt();
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
	if (opt->output_old && opt->n_threads > 1) {
		fprintf(stderr, "[prophyle_index:%s] -v option can be used only with one thread (-t 1)\n", __func__);
		return 1;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "[prophyle_index] Usage:   prophyle_index query [options] <prefix> <in.fq>\n\n");
		fprintf(stderr, "[prophyle_index] Options: -k INT    length of k-mer\n");
		fprintf(stderr, "[prophyle_index]          -u        use klcp for querying\n");
		fprintf(stderr, "[prophyle_index]          -v        output set of chromosomes for every k-mer\n");
		fprintf(stderr, "[prophyle_index]          -p        do not check whether k-mer is on border of two contigs, and show such k-mers in output\n");
		fprintf(stderr, "[prophyle_index]          -l char*  log file name to output statistics\n");
		fprintf(stderr, "[prophyle_index]          -t INT    number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "[prophyle_index] \n");
		return 1;
	}
	if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
		fprintf(stderr, "[prophyle_index:%s] fail to locate the index %s\n", __func__, argv[optind]);
		free(opt);
		return 1;
	}
	query(prefix, argv[optind+1], opt);
	free(opt); free(prefix);
	return 0;
}

int prophyle_index_build(int argc, char *argv[])
{
	int c;
	prophyle_index_opt_t *opt;
	char *prefix;
	opt = prophyle_index_init_opt();
	int sa_intv = 32;
	while ((c = getopt(argc, argv, "si:k:")) >= 0) {
		switch (c) {
		case 'k': opt->kmer_length = atoi(optarg); break;
		case 'i': sa_intv = atoi(optarg); break;
		case 's': opt->construct_sa_parallel = 1; break;
		default: return 1;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "[prophyle_index]\n");
		fprintf(stderr, "[prophyle_index] Usage:   prophyle_index build <prefix>\n\n");
		fprintf(stderr, "[prophyle_index] Options:  -k INT    length of k-mer\n");
		fprintf(stderr, "[prophyle_index]           -s        construct klcp and sa in parallel\n");
		fprintf(stderr, "[prophyle_index]           -i        sampling distance for SA\n");
    fprintf(stderr, "\n");
		return 1;
	}
	if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
		fprintf(stderr, "[prophyle_index:%s] fail to locate the index %s\n", __func__, argv[optind]);
		return 1;
	}
	build_index(prefix, opt, sa_intv);
	free(prefix);
	return 0;
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "[prophyle_index] Program: prophyle_index (alignment of k-mers)\n");
	fprintf(stderr, "[prophyle_index] Usage:   prophyle_index command [options]\n\n");
	fprintf(stderr, "[prophyle_index] Command: build         construct index\n");
	fprintf(stderr, "[prophyle_index] Command: query         query reads against index\n");
	fprintf(stderr, "[prophyle_index] \n");
	return 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	if (argc < 2) return usage();
	if (strcmp(argv[1], "build") == 0) ret = prophyle_index_build(argc - 1, argv + 1);
	else if (strcmp(argv[1], "query") == 0) ret = prophyle_index_query(argc-1, argv+1);
	else return usage();

	return ret;
}
