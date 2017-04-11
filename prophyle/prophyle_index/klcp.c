#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <pthread.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwa_utils.h"
#include "prophyle_query.h"
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"
#include "bwa.h"
#include "bwase.h"
#include "kstring.h"
#include "klcp.h"

uint64_t overall_increase = 0;
int32_t position_of_smallest_zero_bit[MAX_BITARRAY_BLOCK_VALUE + 1];
int32_t position_of_biggest_zero_bit[MAX_BITARRAY_BLOCK_VALUE + 1];

void destroy_klcp(klcp_t* klcp) {
	if (klcp == 0) {
		return;
	}
	destroy_bitarray(klcp->klcp);
	free(klcp);
}

uint64_t decrease_k(klcp_t* klcp, const uint64_t k) {
	int64_t new_k = (int64_t)k;
  new_k = (int64_t)k - 1;
  int stop = 0;
  bitarray_block_t value = klcp->klcp->blocks[new_k / BITS_IN_BLOCK];
	int64_t new_k_res = new_k % BITS_IN_BLOCK;
  value = value >> (BITS_IN_BLOCK - 1 - new_k_res);
  if (value == (1 << (new_k_res + 1)) - 1) {
    new_k -= (new_k_res + 1);
  } else {
    new_k -= position_of_smallest_zero_bit[value];
    stop = 1;
  }
  if (!stop) {
  	while (new_k >= 0) {
  		bitarray_block_t value = klcp->klcp->blocks[new_k / BITS_IN_BLOCK];
      if (value == MAX_BITARRAY_BLOCK_VALUE) {
        new_k -= BITS_IN_BLOCK;
      } else {
        new_k -= position_of_smallest_zero_bit[value];
        break;
      }
    }
  }
  new_k++;
	return (uint64_t)new_k;
}

uint64_t increase_l(klcp_t* klcp, const uint64_t l) {
	int64_t new_l = (int64_t)l;
  new_l = (int64_t)l;
  int stop = 0;
  bitarray_block_t value = klcp->klcp->blocks[new_l / BITS_IN_BLOCK];
  int64_t shift = BITS_IN_BLOCK - new_l % BITS_IN_BLOCK;
  value = value & ((1 << shift) - 1);
  if (value == (1 << shift) - 1) {
    new_l += shift;
  } else {
    new_l += shift - 1 -
      position_of_biggest_zero_bit[(bitarray_block_t)((1 << BITS_IN_BLOCK) - (1 << shift) + value)];
    stop = 1;
  }
  if (!stop) {
	   while (new_l < klcp->seq_len) {
  		bitarray_block_t value = klcp->klcp->blocks[new_l / BITS_IN_BLOCK];
      if (value == MAX_BITARRAY_BLOCK_VALUE) {
        new_l += BITS_IN_BLOCK;
      } else {
        new_l += BITS_IN_BLOCK - 1 - position_of_biggest_zero_bit[value];//find_biggest_zero_index(value);
        stop = 1;
        break;
      }
  	}
  }
  if (new_l > klcp->seq_len) {
    new_l = klcp->seq_len;
  }
	return (uint64_t)new_l;
}

void construct_klcp_recursion(const bwt_t* bwt, bwtint_t k, bwtint_t l, int i, int kmer_length, klcp_t* klcp) {
	if (k > l) {
		return;
	}
	if (k == l) {
		return;
	}
	if (i == kmer_length - 1) {
		uint64_t t;
		for(t = k; t < l; ++t) {
			//klcp->klcp[t] = 1;
			add_to_bitarray(klcp->klcp, t);
		}
		return;
	}
	ubyte_t c = 0;
	bwtint_t new_k = 0;
	bwtint_t new_l = 0;
	bwt_2occ(bwt, k - 1, l, c, &new_k, &new_l);
	new_k = bwt->L2[c] + new_k + 1;
	new_l = bwt->L2[c] + new_l;
	construct_klcp_recursion(bwt, new_k, new_l, i + 1, kmer_length, klcp);
	c++;
	bwt_2occ(bwt, k - 1, l, c, &new_k, &new_l);
	new_k = bwt->L2[c] + new_k + 1;
	new_l = bwt->L2[c] + new_l;
	construct_klcp_recursion(bwt, new_k, new_l, i + 1, kmer_length, klcp);
	c++;
	bwt_2occ(bwt, k - 1, l, c, &new_k, &new_l);
	new_k = bwt->L2[c] + new_k + 1;
	new_l = bwt->L2[c] + new_l;
	construct_klcp_recursion(bwt, new_k, new_l, i + 1, kmer_length, klcp);
	c++;
	bwt_2occ(bwt, k - 1, l, c, &new_k, &new_l);
	new_k = bwt->L2[c] + new_k + 1;
	new_l = bwt->L2[c] + new_l;
	construct_klcp_recursion(bwt, new_k, new_l, i + 1, kmer_length, klcp);
}

void klcp_dump(const char *fn, const klcp_t* klcp)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	err_fwrite(&klcp->seq_len, sizeof(uint64_t), 1, fp);
	err_fwrite(klcp->klcp->blocks, sizeof(bitarray_block_t), klcp->klcp->capacity, fp);
	err_fflush(fp);
	err_fclose(fp);
}

static bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a)
{
	const int bufsize = 0x1000000; // 16M block
	bwtint_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	return offset;
}

int32_t find_smallest_zero_index(bitarray_block_t value) {
  int32_t position = 0;
  while (position < BITS_IN_BLOCK) {
    if (value % 2 != 0) {
      position++;
      value /= 2;
    } else {
      break;
    }
  }
  return position;
}

int32_t find_biggest_zero_index(bitarray_block_t value) {
  int32_t position = BITS_IN_BLOCK - 1;
  while (position >= 0) {
    if ((value & (1 << position)) != 0) {
      position--;
    } else {
      break;
    }
  }
  return position;
}

void klcp_restore(const char *fn, klcp_t* klcp)
{
	FILE *fp;
	fp = xopen(fn, "rb");
	err_fread_noeof(&klcp->seq_len, sizeof(uint64_t), 1, fp);
	klcp->klcp->size = klcp->seq_len;
	klcp->klcp->capacity = (klcp->seq_len + BITS_IN_BLOCK - 1) / BITS_IN_BLOCK;
	klcp->klcp->blocks = (bitarray_block_t*)calloc(klcp->klcp->capacity, sizeof(bitarray_block_t));
	fread_fix(fp, sizeof(bitarray_block_t) * klcp->klcp->capacity, klcp->klcp->blocks);
	err_fclose(fp);
	fprintf(stderr, "klcp was read\n");
  uint64_t i;
  for(i = 0; i <= MAX_BITARRAY_BLOCK_VALUE; ++i) {
    position_of_smallest_zero_bit[i] = find_smallest_zero_index((bitarray_block_t)i);
		position_of_biggest_zero_bit[i] = find_biggest_zero_index((bitarray_block_t)i);
  }
}

klcp_t* construct_klcp(const bwt_t *bwt, const int kmer_length) {
	double t_real;
	t_real = realtime();
	uint64_t n = bwt->seq_len;
  klcp_t* klcp = malloc(sizeof(klcp_t));
  klcp->seq_len = n;
	klcp->klcp = create_bitarray(n);
	uint64_t i;
	for(i = 0; i < klcp->klcp->capacity; ++i) {
		klcp->klcp->blocks[i] = 0;
	}
	construct_klcp_recursion(bwt, (bwtint_t)0, (bwtint_t)n, 0, kmer_length, klcp);
	fprintf(stdout, "\n[%s]  time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	return klcp;
}

typedef struct {
	klcp_t* klcp;
	bwt_t* bwt;
	int kmer_length;
	const char* prefix;
	int sa_intv;
} klcp_data_t;

void* construct_klcp_parallel(void* data) {
	klcp_data_t* klcp_data = (klcp_data_t*)data;
	klcp_data->klcp = construct_klcp(klcp_data->bwt, klcp_data->kmer_length);
	return 0;
}

void* construct_sa_parallel(void* data) {
	klcp_data_t* klcp_data = (klcp_data_t*)data;
	bwt_cal_sa(klcp_data->bwt, klcp_data->sa_intv);
	fprintf(stderr, "sa calculated\n");
	char* fn = malloc((strlen(klcp_data->prefix) + 10) * sizeof(char));
  strcpy(fn, klcp_data->prefix);
  strcat(fn, ".sa");
	bwt_dump_sa(fn, klcp_data->bwt);
	fprintf(stderr, "sa dumped");
	return 0;
}

void exk_index_core(const char *prefix, const exk_opt_t *opt, int sa_intv) {
	bwt_t *bwt;
	{ // load BWT
		if ((bwt = bwa_idx_load_bwt_without_sa(prefix)) == 0) {
			fprintf(stderr, "Couldn't load idx from %s\n", prefix);
			return;
		}
	}

	klcp_t* klcp;
	if (opt->construct_sa_parallel) {
		klcp_data_t* klcp_data = malloc(sizeof(klcp_data_t));
		klcp_data->bwt = bwt;
		klcp_data->kmer_length = opt->kmer_length;
		klcp_data->prefix = prefix;
		klcp_data->sa_intv = sa_intv;
		pthread_t tid[2];
		int status_klcp = pthread_create(&tid[0], NULL, construct_klcp_parallel, (void*)klcp_data);
		int status_sa = pthread_create(&tid[1], NULL, construct_sa_parallel, (void*)klcp_data);
		xassert(!status_klcp, "error while creating thread for klcp parallel construction, try construction separate from sa\n");
		xassert(!status_sa, "error while creating thread for sa parallel construction, try construction separate from klcp\n");
	  fprintf(stderr, "parallel construction of klcp and sa started\n");
		int status_addr_klcp = pthread_join(tid[0], (void**)&status_addr_klcp);
		int status_addr_sa = pthread_join(tid[1], (void**)&status_addr_sa);
		xassert(!status_addr_klcp, "error while klcp parallel construction, try construction separate from sa\n");
		xassert(!status_addr_sa, "error sa parallel construction, try construction separate from klcp\n");
		klcp = klcp_data->klcp;
	} else {
		klcp = construct_klcp(bwt, opt->kmer_length);
	}
  fprintf(stdout, "klcp constructed\n");
  char* fn = malloc((strlen(prefix) + 10) * sizeof(char));
  strcpy(fn, prefix);
  strcat(fn, ".");
  char* kmer_length_str = malloc(5 * sizeof(char));
  sprintf(kmer_length_str, "%d", opt->kmer_length);
  strcat(fn, kmer_length_str);
  strcat(fn, ".klcp");
	klcp_dump(fn, klcp);
  fprintf(stdout, "klcp dumped\n");
	// destroy
	if (opt->construct_sa_parallel) {
		bwt_destroy(bwt);
	} else {
		bwt_destroy_without_sa(bwt);
	}
}
