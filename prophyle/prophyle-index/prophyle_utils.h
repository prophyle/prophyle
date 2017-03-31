#ifndef PROPHYLE_UTILS_H
#define PROPHYLE_UTILS_H

#include <stdint.h>

typedef struct {
	int mode; // bit 24-31 are the barcode length
	int n_threads;
	int trim_qual;

	int use_klcp;
	int kmer_length;
	int output;
	int output_old;
	int output_read_qual;
	int skip_after_fail;
	int skip_positions_on_border;

	int need_log;
	char* log_file_name;

	int construct_sa_parallel;
} exk_opt_t;

#endif //
