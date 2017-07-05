#pragma once

#include "tree_index.h"
#include <vector>
#include <iostream>
#include <set>
#include <sstream>

enum class AssignmentOutputFormat {
  Sam = 0,
  Kraken,
  Count
};

enum class Measure {
  H1 = 0,
  C1,
  Count
};

class ReadProcessor {
public:
  ReadProcessor(const TreeIndex& tree, size_t k, bool simulate_lca = false, bool annotate = false,
      bool tie_lca = false, bool not_translate_blocks = false);

  void process_krakline(const std::string& krakline, AssignmentOutputFormat format, Measure criteria);
  void print_sam_header(std::ostream& out) const;

private:
  static constexpr size_t kFakeContigLength = 42424242;

  void load_krakline(const std::string& krakline);
  void filter_assignments();
  void print_assignments(AssignmentOutputFormat format, Measure criteria);
  void print_sam_line(int32_t node_id, const std::string& suffix = "", std::ostream& out = std::cout);
  void print_kraken_line(int32_t node_id, std::ostream& out = std::cout);
  void fill_masks_from_kmer_blocks();
  void set_masks(const std::vector<int32_t>& node_ids, size_t block_position, size_t block_length);
  void print_masks() const;
  void clear();
  void clear_masks();
  void propagate_matching_kmers();
  void copy_masks(int32_t node_from, int32_t node_to);
  void fill_cigar(const std::vector<int16_t>& mask, size_t mask_size, std::string& cigar);

  size_t hit_mask_size() const {
    return read_length_ - k_ + 1;
  }

  size_t coverage_mask_size() const {
    return read_length_;
  }

  const TreeIndex& tree_;
  size_t k_;
  const bool simulate_lca_;
  const bool annotate_;
  const bool tie_lca_;
  const bool not_translate_blocks_;

  std::string read_name_;
  size_t read_length_;
  std::string krakmers_;
  std::string read_;
  std::string qualities_;

  std::vector<std::vector<int16_t>> hit_masks_;
  std::vector<std::vector<int16_t>> coverage_masks_;
  std::set<int32_t> matching_nodes_;

  std::vector<int32_t> best_hit_nodes_;
  int32_t best_hit_;
  std::vector<int32_t> best_coverage_nodes_;
  int32_t best_coverage_;

  std::string hit_cigar_;
  std::string coverage_cigar_;
};
