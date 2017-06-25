#include "tree_index.h"
#include "knhx.h"
#include <sstream>
#include <fstream>

TreeIndex::TreeIndex(const std::string& tree_filename) {
  std::ifstream t(tree_filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  const char *s = buffer.str().c_str();
  int error;
  first_node_ = kn_parse(s, &nodes_count_, &error);
  bool root_found = false;
  for (int32_t i = 0; i < nodes_count_; ++i) {
    knhx1_t* possible_root = first_node_ + i;
    if (possible_root->parent == -1) {
      root_ = possible_root;
      root_found = true;
      break;
    }
  }
  if (!root_found) {
    std::cerr << "root not found" << std::endl;
    exit(1);
  }
  for (int32_t i = 0; i < nodes_count_; ++i) {
    const knhx1_t* node = node_by_id(i);
    id_by_name_.emplace(std::pair<std::string, int32_t>(node->name, i));
  }
  upper_nodes_.resize(nodes_count_);
  fill_upper_nodes(root_);
}

void TreeIndex::fill_upper_nodes(const knhx1_t* node) {
  int32_t node_id = id_by_name(node->name);
  // std::cerr << "processing node: " << node_id << " " << node->name << " " << node->n << std::endl;
  for (int32_t i = 0; i < node->n; ++i) {
    auto child = node_by_id(node->child[i]);
    upper_nodes_[node->child[i]] = upper_nodes_[node_id];
    upper_nodes_[node->child[i]].push_back(node_id);
    fill_upper_nodes(child);
  }
}

std::ostream& operator<<(std::ostream& sout, const TreeIndex& tree) {
  kstring_t str;
  for (int32_t i = 0; i < tree.nodes_count(); ++i) {
    const knhx1_t *p = tree.node_by_id(i);
    printf("[%d] %s\t%d\t%d\t%g", i, p->name, p->parent, p->n, p->d);
    for (int32_t j = 0; j < p->n; ++j) {
      printf("\t%d", p->child[j]);
    }
    putchar('\n');
  }
  str.l = str.m = 0; str.s = 0;
  kn_format(tree.first_node(), tree.nodes_count() - 1, &str);
  puts(str.s);
  free(str.s);
  return sout;
}
