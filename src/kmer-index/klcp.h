#ifndef KLCP_H
#define KLCP_H

#include <stdint.h>
#include "bwt.h"
#include "exkutils.h"
#include "bitarray.h"

typedef struct {
  int seq_len;
  bitarray_t* klcp;
  uint64_t* prev;
  uint64_t* next;
} klcp_t;

klcp_t* construct_klcp(const bwt_t *bwt, const int kmer_length);
void klcp_restore(const char *fn, klcp_t* klcp);
void exk_index_core(const char *prefix, const char *fn_fa, const exk_opt_t *opt);
uint64_t decrease_k(klcp_t* klcp, const uint64_t k);
uint64_t increase_l(klcp_t* klcp, const uint64_t l);

#endif //KLCP_H
