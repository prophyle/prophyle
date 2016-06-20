#include<stdio.h>
#include<stdlib.h>
#include <stdint.h>

typedef struct {
  char* values;
  uint64_t size;
  uint64_t capacity;
} bitarray_t;

bitarray_t* create_bitarray(uint64_t n);

void add_to_bitarray(bitarray_t* array, uint64_t value);

void delete_from_bitarray(bitarray_t* array, uint64_t value);

int is_member(bitarray_t* array, uint64_t value);
