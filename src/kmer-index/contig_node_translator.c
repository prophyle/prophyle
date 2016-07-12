#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "contig_node_translator.h"

int contig_to_node[200000000];
char* node_names[100000];
int nodes_count = 0;
int contigs_count = 0;

int get_node_from_contig(int contig) {
  if (contig < 0 || contig >= contigs_count) {
    fprintf(stderr, "contig %d is outside of range [%d, %d]\n", contig, 0, contigs_count - 1);
  }
  return contig_to_node[contig];
}

void add_contig(const char* contig, int contig_number) {
  contigs_count++;
  char* ch = strchr(contig, '_');
  int index = ch - contig;
  char* node_name = malloc((index + 1) * sizeof(char));
  memcpy(node_name, contig, index);
  node_name[index] = '\0';
  // fprintf(stderr, "node_name = %s\n", node_name);
  // if (nodes_count > 0) {
  //   fprintf(stderr, "prev node_name = %s\n", node_names[nodes_count - 1]);
  //   fprintf(stderr, "equal = %d\n", strcmp(node_name, node_names[nodes_count - 1]));
  // }
  if (nodes_count == 0 || strcmp(node_name, node_names[nodes_count - 1])) {
    node_names[nodes_count] = node_name;
    contig_to_node[contig_number] = nodes_count;
    nodes_count++;
  } else {
    contig_to_node[contig_number] = nodes_count - 1;
  }
  // fprintf(stderr, "contig = %d, node = %d\n", contig_number, contig_to_node[contig_number]);
}
