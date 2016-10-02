#include<stdio.h>
#include<stdlib.h>
#include <stdint.h>

#define bitarray_value_t uint16_t
#define BITS_IN_VALUE 16
#define MAX_BITARRAY_VALUE (1 << BITS_IN_VALUE) - 1

typedef struct {
  bitarray_value_t* values;
  uint64_t size;
  uint64_t capacity;
} bitarray_t;

void destroy_bitarray(bitarray_t* array);
bitarray_t* create_bitarray(uint64_t n);
void add_to_bitarray(bitarray_t* array, uint64_t value);
void delete_from_bitarray(bitarray_t* array, uint64_t value);
int is_member(bitarray_t* array, uint64_t value);
