#include<stdio.h>
#include<stdlib.h>

typedef struct {
  char* values;
  int size;
  int capacity;
} bitarray_t;

bitarray_t* create_bitarray(int n);

void add_to_bitarray(bitarray_t* array, int value);

void delete_from_bitarray(bitarray_t* array, int value);

int is_member(bitarray_t* array, int value);
