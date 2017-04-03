#ifndef KLCP_H
#define KLCP_H

#include <stdint.h>
#include "bwt.h"
#include "prophyle_utils.h"
#include "bitarray.h"

typedef struct {
  uint64_t seq_len;
  bitarray_t* klcp;
} klcp_t;

void destroy_klcp(klcp_t* klcp);
klcp_t* construct_klcp(const bwt_t *bwt, const int kmer_length);
void klcp_restore(const char *fn, klcp_t* klcp);
void exk_index_core(const char *prefix, const exk_opt_t *opt, int sa_intv);
uint64_t decrease_k(klcp_t* klcp, const uint64_t k);
uint64_t increase_l(klcp_t* klcp, const uint64_t l);

#endif //KLCP_H
