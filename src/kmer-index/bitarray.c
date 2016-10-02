#include<stdio.h>
#include<stdlib.h>
#include "bitarray.h"

void destroy_bitarray(bitarray_t* array) {
  if (array == 0) {
    return;
  }
  free(array->values);
  free(array);
}

bitarray_t* create_bitarray(uint64_t n)
{
  bitarray_t* array = malloc(sizeof(bitarray_t));
  array->size = n;
  array->capacity = (n + BITS_IN_VALUE - 1) / BITS_IN_VALUE;
  array->values = calloc(array->capacity, sizeof(bitarray_value_t));
  return array;
}

void add_to_bitarray(bitarray_t* array, uint64_t value)
{
  array->values[value / BITS_IN_VALUE] =
    array->values[value / BITS_IN_VALUE] | ( 1 << (BITS_IN_VALUE - 1 - value % BITS_IN_VALUE));
}

void delete_from_bitarray(bitarray_t* array, uint64_t value)
{
  array->values[value / BITS_IN_VALUE] =
    array->values[value / BITS_IN_VALUE] & ~(1 << (BITS_IN_VALUE - 1 - value % BITS_IN_VALUE));
}

int is_member(bitarray_t* array, uint64_t value)
{
  return (array->values[value / BITS_IN_VALUE] & (1 << (BITS_IN_VALUE - 1 - value % BITS_IN_VALUE)));
}
