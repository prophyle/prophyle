#ifndef PROPHYLE_QUERY_H
#define PROPHYLE_QUERY_H

#include <stdint.h>
#include "bwt.h"
#include "bwtaln.h"
#include "bwa.h"
#include "klcp.h"
#include "prophyle_utils.h"

typedef struct {
	uint64_t position;
	int strand;
	int rid;
	int node;
} bwt_position_t;

void prophyle_index_query_core(const char* prefix, const char* fn_fa, const prophyle_index_opt_t* opt);

#endif //PROPHYLE_QUERY_H
