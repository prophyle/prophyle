#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "contig_node_translator.h"

int contig_to_node[200000000];
char* node_names[100000];
int nodes_count = 0;

void add_contig(const char* contig, int contig_number) {
  char* ch = strchr(contig, '_');
  int index = ch - contig;
  char* node_name = malloc((index + 1) * sizeof(char));
  memcpy(node_name, contig, index);
  node_name[index] = '\0';
  fprintf(stderr, "node_name = %s\n", node_name);
  if (nodes_count > 0) {
    fprintf(stderr, "prev node_name = %s\n", node_names[nodes_count - 1]);
    fprintf(stderr, "equal = %d\n", strcmp(node_name, node_names[nodes_count - 1]));
  }
  if (nodes_count == 0 || strcmp(node_name, node_names[nodes_count - 1])) {
    node_names[nodes_count] = node_name;
    contig_to_node[contig_number] = nodes_count;
    nodes_count++;
  } else {
    contig_to_node[contig_number] = nodes_count - 1;
  }
  fprintf(stderr, "contig = %d, node = %d\n", contig_number, contig_to_node[contig_number]);
}
