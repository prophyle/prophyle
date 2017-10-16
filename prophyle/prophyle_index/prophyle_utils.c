#include <stdlib.h>
#include "prophyle_utils.h"

prophyle_index_opt_t* prophyle_index_init_opt()
{
	prophyle_index_opt_t* o;
	o = (prophyle_index_opt_t*)calloc(1, sizeof(prophyle_index_opt_t));
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->n_threads = 1;
	o->trim_qual = 0;
	o->kmer_length = 14;
	o->use_klcp = 0;
	o->output = 1;
	o->output_read_qual = 0;
	o->output_old = 0;
	o->skip_positions_on_border = 1;
	o->construct_sa_parallel = 0;
	o->need_log = 0;
	o->log_file_name = NULL;
	return o;
}
