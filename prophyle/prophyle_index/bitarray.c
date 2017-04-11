#include <stdio.h>
#include <stdlib.h>
#include "bitarray.h"

void destroy_bitarray(bitarray_t* array) {
  if (array == 0) {
    return;
  }
  free(array->blocks);
  free(array);
}

bitarray_t* create_bitarray(uint64_t n)
{
  bitarray_t* array = malloc(sizeof(bitarray_t));
  array->size = n;
  array->capacity = (n + BITS_IN_BLOCK - 1) / BITS_IN_BLOCK;
	array->blocks = calloc(array->capacity, sizeof(bitarray_block_t));
  return array;
}

void add_to_bitarray(bitarray_t* array, uint64_t value)
{
	array->blocks[value / BITS_IN_BLOCK] =
    array->blocks[value / BITS_IN_BLOCK] | ( 1 << (BITS_IN_BLOCK - 1 - value % BITS_IN_BLOCK));
}

void delete_from_bitarray(bitarray_t* array, uint64_t value)
{
	array->blocks[value / BITS_IN_BLOCK] =
    array->blocks[value / BITS_IN_BLOCK] & ~(1 << (BITS_IN_BLOCK - 1 - value % BITS_IN_BLOCK));
}
