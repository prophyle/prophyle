#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "utils.h"
#include "klcp.h"

int32_t position_of_smallest_zero_bit[MAX_BITARRAY_BLOCK_VALUE + 1];
int32_t position_of_biggest_zero_bit[MAX_BITARRAY_BLOCK_VALUE + 1];

void destroy_klcp(klcp_t* klcp) {
	if (klcp == 0) {
		return;
	}
	destroy_bitarray(klcp->klcp);
	free(klcp);
}

uint64_t decrease_sa_position(const klcp_t* klcp, uint64_t k) {
	int64_t new_position = (int64_t)k - 1;
  int new_position_found = 0;
  bitarray_block_t value = klcp->klcp->blocks[new_position / BITS_IN_BLOCK];
	int64_t new_position_res = new_position % BITS_IN_BLOCK;
  value = value >> (BITS_IN_BLOCK - 1 - new_position_res);
  if (value == (1 << (new_position_res + 1)) - 1) {
    new_position -= (new_position_res + 1);
  } else {
    new_position -= position_of_smallest_zero_bit[value];
    new_position_found = 1;
  }
  if (!new_position_found) {
  	while (new_position >= 0) {
  		bitarray_block_t value = klcp->klcp->blocks[new_position / BITS_IN_BLOCK];
      if (value == MAX_BITARRAY_BLOCK_VALUE) {
        new_position -= BITS_IN_BLOCK;
      } else {
        new_position -= position_of_smallest_zero_bit[value];
        break;
      }
    }
  }
  new_position++;
	return (uint64_t)new_position;
}

uint64_t increase_sa_position(const klcp_t* klcp, uint64_t l) {
	int64_t new_position = (int64_t)l;
  int new_position_found = 0;
  bitarray_block_t value = klcp->klcp->blocks[new_position / BITS_IN_BLOCK];
  int64_t shift = BITS_IN_BLOCK - new_position % BITS_IN_BLOCK;
  value = value & ((1 << shift) - 1);
  if (value == (1 << shift) - 1) {
    new_position += shift;
  } else {
    new_position += shift - 1 -
      position_of_biggest_zero_bit[(bitarray_block_t)((1 << BITS_IN_BLOCK) - (1 << shift) + value)];
    new_position_found = 1;
  }
  if (!new_position_found) {
	   while (new_position < klcp->seq_len) {
  		bitarray_block_t value = klcp->klcp->blocks[new_position / BITS_IN_BLOCK];
      if (value == MAX_BITARRAY_BLOCK_VALUE) {
        new_position += BITS_IN_BLOCK;
      } else {
        new_position += BITS_IN_BLOCK - 1 - position_of_biggest_zero_bit[value];
        new_position_found = 1;
        break;
      }
  	}
  }
  if (new_position > klcp->seq_len) {
    new_position = klcp->seq_len;
  }
	return (uint64_t)new_position;
}

void construct_klcp_recursion(const bwt_t* bwt, bwtint_t k, bwtint_t l, int tree_depth, int kmer_length, klcp_t* klcp) {
	if (k > l) {
		return;
	}
	if (k == l) {
		return;
	}
	if (tree_depth == kmer_length - 1) {
		uint64_t sa_position;
		for(sa_position = k; sa_position < l; ++sa_position) {
			add_to_bitarray(klcp->klcp, sa_position);
		}
		return;
	}
	ubyte_t c = 0;
	bwtint_t new_k = 0;
	bwtint_t new_l = 0;
	bwt_2occ(bwt, k - 1, l, c, &new_k, &new_l);
	construct_klcp_recursion(bwt, bwt->L2[c] + new_k + 1, bwt->L2[c] + new_l, tree_depth + 1, kmer_length, klcp);
	c++;
	bwt_2occ(bwt, k - 1, l, c, &new_k, &new_l);
	construct_klcp_recursion(bwt, bwt->L2[c] + new_k + 1, bwt->L2[c] + new_l, tree_depth + 1, kmer_length, klcp);
	c++;
	bwt_2occ(bwt, k - 1, l, c, &new_k, &new_l);
	construct_klcp_recursion(bwt, bwt->L2[c] + new_k + 1, bwt->L2[c] + new_l, tree_depth + 1, kmer_length, klcp);
	c++;
	bwt_2occ(bwt, k - 1, l, c, &new_k, &new_l);
	construct_klcp_recursion(bwt, bwt->L2[c] + new_k + 1, bwt->L2[c] + new_l, tree_depth + 1, kmer_length, klcp);
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
    if ((value & (1 << position)) != 0) {
      position++;
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
	fprintf(stderr, "[prophyle_index:%s]  time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	return klcp;
}
