/*
  KLCP data structure block-based implementation.
  Author: Kamil Salikhov <salikhov.kamil@gmail.com>
  Licence: MIT
*/

#ifndef KLCP_H
#define KLCP_H

#include "bwt.h"
#include "bitarray.h"

typedef struct {
  uint64_t seq_len;
  bitarray_t* klcp;
} klcp_t;

void destroy_klcp(klcp_t* klcp);
void klcp_dump(const char *fn, const klcp_t* klcp);
klcp_t* construct_klcp(const bwt_t *bwt, const int kmer_length);
void klcp_restore(const char *fn, klcp_t* klcp);
uint64_t decrease_sa_position(const klcp_t* klcp, uint64_t position);
uint64_t increase_sa_position(const klcp_t* klcp, uint64_t position);

#endif //KLCP_H
