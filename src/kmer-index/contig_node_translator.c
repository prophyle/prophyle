#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "contig_node_translator.h"

static int contig_to_node[200000000];
static char* node_names[100000];
static int nodes_count = 0;
static int contigs_count = 0;

int get_node_from_contig(int contig) {
  if (contig < 0 || contig >= contigs_count) {
    fprintf(stderr, "contig %d is outside of range [%d, %d]\n", contig, 0, contigs_count - 1);
  }
  return contig_to_node[contig];
}

char* get_node_name(const int node) {
  return node_names[node];
}

void add_contig(const char* contig, int contig_number) {
  contigs_count++;
  char* ch = strchr(contig, '_');
  int index = 0;
  if (ch == NULL) {
    index = strlen(contig);
  } else {
    index = ch - contig;
  }
  char* node_name = malloc((index + 1) * sizeof(char));
  memcpy(node_name, contig, index);
  node_name[index] = '\0';
  if (nodes_count == 0 || strcmp(node_name, node_names[nodes_count - 1])) {
    node_names[nodes_count] = node_name;
    contig_to_node[contig_number] = nodes_count;
    nodes_count++;
  } else {
    contig_to_node[contig_number] = nodes_count - 1;
    free(node_name);
  }
  //fprintf(stderr, "contig = %d, nodes_count = %d\n", contig_number, nodes_count);
}
