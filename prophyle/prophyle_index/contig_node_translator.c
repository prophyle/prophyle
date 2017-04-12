#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "contig_node_translator.h"
#include "utils.h"

#define MAX_NODES_COUNT 100000
#define MAX_CONTIGS_COUNT 200000000

static int contig_to_node[MAX_CONTIGS_COUNT];
static char* node_names[MAX_NODES_COUNT];
static int node_name_lengths[MAX_NODES_COUNT];
static int nodes_count = 0;
static int contigs_count = 0;

int get_node_from_contig(int contig) {
  if (contig < 0 || contig >= contigs_count) {
    fprintf(stderr, "[prophyle_index:%s] contig %d is outside of range [%d, %d]\n",
      __func__, contig, 0, contigs_count - 1);
  }
  return contig_to_node[contig];
}

char* get_node_name(int node) {
  return node_names[node];
}

int get_node_name_length(int node) {
  return node_name_lengths[node];
}

void add_contig(char* contig, int contig_number) {
  xassert(contigs_count < MAX_CONTIGS_COUNT,
    "[prophyle_index] there are more than MAX_CONTIGS_COUNT contigs, try to increase MAX_CONTIGS_COUNT in contig_node_translator.c\n");
  contigs_count++;
  const char* ch = strchr(contig, '@');
  int index = 0;
  if (ch == NULL) {
    index = strlen(contig) - 1;
  } else {
    index = ch - contig;
  }
  contig[index] = '\0';
  if (nodes_count == 0 || strcmp(contig, node_names[nodes_count - 1])) {
    char* node_name = malloc((index + 1) * sizeof(char));
    memcpy(node_name, contig, index);
    node_name[index] = '\0';
    xassert(nodes_count < MAX_NODES_COUNT,
      "[prophyle_index] there are more than MAX_NODES_COUNT nodes, try to increase MAX_NODES_COUNT in contig_node_translator.c\n");
    node_names[nodes_count] = node_name;
    node_name_lengths[nodes_count] = strlen(node_name);
    contig_to_node[contig_number] = nodes_count;
    nodes_count++;
  } else {
    contig_to_node[contig_number] = nodes_count - 1;
  }
}
