#include<stdio.h>
#include<stdlib.h>
#include "bitarray.h"

bitarray_t* create_bitarray(uint64_t n)
{
  bitarray_t* array = malloc(sizeof(bitarray_t));
  array->size = n;
  array->capacity = (n + 7) / 8;
	array->values = calloc(array->capacity, sizeof(char));
  return array;
}

void add_to_bitarray(bitarray_t* array, uint64_t value)
{
	array->values[value / 8] = array->values[value / 8] | ( 1 << (value % 8));
}

void delete_from_bitarray(bitarray_t* array, uint64_t value)
{
	array->values[value / 8] = array->values[value / 8] & ~(1 << (value % 8));
}

int is_member(bitarray_t* array, uint64_t value)
{
	return (array->values[value / 8] & (1 << (value % 8)));
}
