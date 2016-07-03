#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtexk.h"
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"
#include "bwa.h"
#include "bwase.h"
#include "kstring.h"
#include "klcp.h"

const uint32_t SAMPLING_DISTANCE = 128;

uint64_t overall_increase = 0;

void destroy_klcp(klcp_t* klcp) {
	if (klcp == 0) {
		return;
	}
	destroy_bitarray(klcp->klcp);
	free(klcp->next);
	free(klcp->prev);
	free(klcp);
}

uint64_t prev_size(const klcp_t* klcp) {
	return (klcp->seq_len + SAMPLING_DISTANCE - 1) / SAMPLING_DISTANCE;
}

uint64_t next_size(const klcp_t* klcp) {
	return (klcp->seq_len + 1 + SAMPLING_DISTANCE - 1) / SAMPLING_DISTANCE;
}

uint64_t decrease_k(klcp_t* klcp, const uint64_t k) {
	uint64_t new_k = k;
	uint64_t size = prev_size(klcp);
	while (new_k >= 1 && is_member(klcp->klcp, new_k - 1)) {
		overall_increase++;
		new_k--;
		if (new_k % SAMPLING_DISTANCE == 0 && new_k < size / SAMPLING_DISTANCE) {
			new_k = klcp->prev[new_k / SAMPLING_DISTANCE];
			break;
		}
	}
	return new_k;
}

uint64_t increase_l(klcp_t* klcp, const uint64_t l) {
	uint64_t new_l = l;
	uint64_t size = next_size(klcp);
	while (new_l < klcp->seq_len && is_member(klcp->klcp, new_l)) {
		overall_increase++;
		new_l++;
		if (new_l % SAMPLING_DISTANCE == 0 && new_l < size / SAMPLING_DISTANCE) {
			new_l = klcp->next[new_l / SAMPLING_DISTANCE];
			break;
		}
	}
	return new_l;
}

void construct_klcp_recursion(const bwt_t* bwt, bwtint_t k, bwtint_t l, int i, int kmer_length, klcp_t* klcp) {
	if (k > l) {
		//fprintf(stderr, "exit on i = %d\n", i);
		return;
	}
	if (k == l) {
		//fprintf(stderr, "exit on i = %d, k = l = %d\n", i, k);
		return;
	}
	if (i == kmer_length - 1) {
		//fprintf(stderr, "New entry found! k = %d, l = %d\n", k, l);
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
	err_fwrite(klcp->klcp->values, sizeof(char), klcp->klcp->capacity, fp);
	err_fwrite(klcp->prev, sizeof(uint64_t), prev_size(klcp), fp);
	err_fwrite(klcp->next, sizeof(uint64_t), next_size(klcp), fp);
	err_fflush(fp);
	err_fclose(fp);
}

static bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a)
{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
	const int bufsize = 0x1000000; // 16M block
	bwtint_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	return offset;
}

void construct_aux_arrays(klcp_t* klcp) {
	int take_current = 1;
	uint64_t prev_one = -1;
	klcp->prev = malloc(prev_size(klcp) * sizeof(uint64_t));
	int64_t i;
	for(i = 0; i < klcp->seq_len; ++i) {
		if (take_current) {
			if (i % SAMPLING_DISTANCE == 0) {
				klcp->prev[i / SAMPLING_DISTANCE] = i;
			}
		} else {
			if (i % SAMPLING_DISTANCE == 0) {
				klcp->prev[i / SAMPLING_DISTANCE] = prev_one;
			}
		}
		if (is_member(klcp->klcp, i)) {
			if (take_current) {
				take_current = 0;
				prev_one = i;
			}
		} else {
			take_current = 1;
		}
	}
	take_current = 1;
	uint64_t next_zero = -1;
	klcp->next = malloc(next_size(klcp) * sizeof(uint64_t));
	for(i = klcp->seq_len; i >= 0; --i) {
		if (i < klcp->seq_len && is_member(klcp->klcp, i)) {
			if (take_current) {
				next_zero = i + 1;
				take_current = 0;
				if (i % SAMPLING_DISTANCE == 0) {
					klcp->next[i / SAMPLING_DISTANCE] = next_zero;
				}
			} else {
				if (i % SAMPLING_DISTANCE == 0) {
					klcp->next[i / SAMPLING_DISTANCE] = next_zero;
				}
			}
		} else {
			if (i % SAMPLING_DISTANCE == 0) {
				klcp->next[i / SAMPLING_DISTANCE] = i;
			}
			take_current = 1;
		}
	}
}

void klcp_restore(const char *fn, klcp_t* klcp)
{
	FILE *fp;
	fp = xopen(fn, "rb");
	err_fread_noeof(&klcp->seq_len, sizeof(uint64_t), 1, fp);
	klcp->klcp->size = klcp->seq_len;
	klcp->klcp->capacity = (klcp->seq_len + 7) / 8;
	klcp->klcp->values = (char*)calloc(klcp->klcp->capacity, sizeof(char));
	klcp->prev = (uint64_t*)calloc(prev_size(klcp), sizeof(uint64_t));
	klcp->next = (uint64_t*)calloc(next_size(klcp), sizeof(uint64_t));
	fread_fix(fp, sizeof(int8_t) * klcp->klcp->capacity, klcp->klcp->values);
	fread_fix(fp, sizeof(uint64_t) * prev_size(klcp), klcp->prev);
	fread_fix(fp, sizeof(uint64_t) * next_size(klcp), klcp->next);
	err_fclose(fp);
	fprintf(stderr, "klcp was read\n");
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
		klcp->klcp->values[i] = 0;
	}
	construct_klcp_recursion(bwt, (bwtint_t)0, (bwtint_t)n, 0, kmer_length, klcp);
	construct_aux_arrays(klcp);
	// fprintf(stderr, "DIRECT\n");
	// construct_klcp_direct(bwt, kmer_length, klcp);

	fprintf(stdout, "\n[%s]  time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	return klcp;
}

void exk_index_core(const char *prefix, const char *fn_fa, const exk_opt_t *opt) {
	bwt_t *bwt;
	// initialization
	{ // load BWT
		//fprintf(stderr, "%s\n", prefix);
		if ((bwt = bwa_idx_load_bwt(prefix)) == 0) {
			fprintf(stderr, "Couldn't load idx from %s\n", prefix);
			return;
		}
	}

	klcp_t* klcp = construct_klcp(bwt, opt->kmer_length);
  fprintf(stdout, "klcp constructed\n");
  char* fn = malloc((strlen(prefix) + 10) * sizeof(char));
  strcpy(fn, prefix);
  strcat(fn, ".");
  char* kmer_length_str = malloc(5 * sizeof(char));
  sprintf(kmer_length_str, "%d", opt->kmer_length);
  strcat(fn, kmer_length_str);
  strcat(fn, ".bit.klcp");
	klcp_dump(fn, klcp);
  fprintf(stdout, "klcp dumped\n");
	// destroy
	bwt_destroy(bwt);
}

// void construct_klcp_direct(const bwt_t* bwt, int kmer_length, klcp_t* klcp) {
// 	int* ks = malloc(kmer_length * sizeof(int));
// 	int* ls = malloc(kmer_length * sizeof(int));
// 	int* cs = malloc(kmer_length * sizeof(ubyte_t));
// 	for(int i = 0; i < kmer_length; ++i) {
// 		ks[i] = 0;
// 		ls[i] = 0;
// 	}
// 	bwtint_t k = 0, l = bwt->seq_len;
// 	ubyte_t c = 0;
// 	int direction = 1;
// 	int i = 0;
// 	while (i >= 0) {
// 		if (i == 0) {
// 			k = 0;
// 			l = bwt->seq_len;
// 		}
// 		else {
// 			k = ks[i - 1];
// 			l = ls[i - 1];
// 		}
// 		if (k > l) {
// 			direction = -1;
// 			i--;
// 			continue;
// 		}
// 		if (k == l) {
// 			direction = -1;
// 			i--;
// 			continue;
// 		}
// 		if (direction == 1) {
// 			if (i == kmer_length) {
// 				//fprintf(stderr, "New entry found! k = %d, l = %d\n", k, l);
// 				// for(int t = 0; t < kmer_length; ++t) {
// 				// 	fprintf(stderr, "%d %d %d\n", ks[t], ls[t], cs[t]);
// 				// }
// 				for(int t = k; t < l; ++t) {
// 					klcp[t] = 1;
// 				}
// 				direction = -1;
// 				i--;
// 				continue;
// 			}
// 			c = 0;
// 			bwtint_t new_k = 0;
// 			bwtint_t new_l = 0;
// 			bwt_2occ(bwt, k - 1, l, c, &new_k, &new_l);
// 			new_k = bwt->L2[c] + new_k + 1;
// 			new_l = bwt->L2[c] + new_l;
// 			ks[i] = new_k;
// 			ls[i] = new_l;
// 			cs[i] = c;
// 			direction = 1;
// 			i++;
// 			continue;
// 		} else {
// 			if (cs[i] == 3) {
// 				direction = -1;
// 				i--;
// 				continue;
// 			}
// 			cs[i]++;
// 			c = cs[i];
// 			bwtint_t new_k = 0;
// 			bwtint_t new_l = 0;
// 			bwt_2occ(bwt, k - 1, l, c, &new_k, &new_l);
// 			new_k = bwt->L2[c] + new_k + 1;
// 			new_l = bwt->L2[c] + new_l;
// 			ks[i] = new_k;
// 			ls[i] = new_l;
// 			direction = 1;
// 			i++;
// 			continue;
// 		}
// 	}
// }
