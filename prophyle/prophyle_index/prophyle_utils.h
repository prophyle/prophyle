/*
	Structure for prophyle_index options.
	Author: Kamil Salikhov <salikhov.kamil@gmail.com>
	Licence: MIT
*/

#ifndef PROPHYLE_UTILS_H
#define PROPHYLE_UTILS_H

#include <stdint.h>
#include <stdlib.h>
#include "bwtaln.h"

typedef struct {
	int mode;
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
} prophyle_index_opt_t;

prophyle_index_opt_t* prophyle_index_init_opt();

#endif //PROPHYLE_UTILS_H
