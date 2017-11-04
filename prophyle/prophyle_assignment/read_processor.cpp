#include "read_processor.h"
#include "word_splitter.h"
#include "knhx.h"
#include <vector>
#include <sstream>
#include <iostream>
#include <string.h>
#include <experimental/string_view>

namespace {

int32_t ones_count(const std::vector<int16_t> &values, size_t size) {
  int32_t count = 0;
  for (size_t i = 0; i < size; ++i) {
    count += values[i];
  }
  return count;
}

}

constexpr size_t ReadProcessor::kFakeContigLength;

ReadProcessor::ReadProcessor(const TreeIndex& tree, size_t k, bool simulate_lca, bool annotate,
    bool tie_lca, bool not_translate_blocks):
    tree_(tree),
    k_(k),
    simulate_lca_(simulate_lca),
    annotate_(annotate),
    tie_lca_(tie_lca),
    not_translate_blocks_(not_translate_blocks),
    hit_masks_(tree_.nodes_count()),
    coverage_masks_(tree_.nodes_count()) {
  if (simulate_lca_) {
    std::cerr << "simulating lca is not supported yet" << std::endl;
    exit(1);
  }
  if (tie_lca_) {
    std::cerr << "tie lca is not supported yet" << std::endl;
    exit(1);
  }
}

void ReadProcessor::process_read(const std::string& krakline, AssignmentOutputFormat format, Measure criteria) {
  parse_krakline(krakline);
  propagate_matching_kmers();
  filter_assignments();
  print_assignments(format, criteria);
  clear();
}

void ReadProcessor::filter_assignments() {
  best_hit_ = -1;
  for (auto node_id : matching_nodes_) {
    auto& hit_mask = hit_masks_[node_id];
    auto hit = ones_count(hit_mask, hit_mask_size());
    if (hit > best_hit_) {
      best_hit_ = hit;
      best_hit_nodes_.clear();
      best_hit_nodes_.push_back(node_id);
    } else if (hit == best_hit_) {
      best_hit_nodes_.push_back(node_id);
    }
  }
  best_coverage_ = -1;
  for (auto node_id : matching_nodes_) {
    auto& coverage_mask = coverage_masks_[node_id];
    auto coverage = ones_count(coverage_mask, coverage_mask_size());
    if (coverage > best_coverage_) {
      best_coverage_ = coverage;
      best_coverage_nodes_.clear();
      best_coverage_nodes_.push_back(node_id);
    } else if (coverage == best_coverage_) {
      best_coverage_nodes_.push_back(node_id);
    }
  }
}

void ReadProcessor::print_assignments(AssignmentOutputFormat format, Measure criteria) {
  std::vector<int32_t> best_matching_nodes;
  if (criteria == Measure::H1) {
    best_matching_nodes = best_hit_nodes_;
  } else {
    best_matching_nodes = best_coverage_nodes_;
  }
  bool tie_solved = false;
  if (tie_lca_ && best_matching_nodes.size() > 1) {
    std::cerr << "tie_lca is not implemented yet" << std::endl;
    exit(1);
  }
  if (best_matching_nodes.size() > 0) {
    for (auto id : best_matching_nodes) {
      if (format == AssignmentOutputFormat::Sam) {
        if (!tie_solved) {
          fill_cigar(hit_masks_[id], hit_mask_size(), hit_cigar_);
          fill_cigar(coverage_masks_[id], coverage_mask_size(), coverage_cigar_);
        }
        // add annotations
        print_sam_line(id, annotate_ ? tree_.joined_tags(id) : "");
      } else if (format == AssignmentOutputFormat::Kraken) {
        print_kraken_line(id);
      }
    }
  } else {
    if (format == AssignmentOutputFormat::Sam) {
      print_sam_line(-1);
    } else if (format == AssignmentOutputFormat::Kraken) {
      print_kraken_line(-1);
    }
  }
}

void ReadProcessor::print_sam_header(std::ostream& out) const {
  out << "@HD\tVN:1.5\tSO:unsorted" << std::endl;
  for (int32_t id = 0; id < tree_.nodes_count(); ++id) {
    const auto& tags = tree_.tags(id);
    auto node = tree_.node_by_id(id);
    auto gi = tags.find("gi");
    auto fastapath = tags.find("fastapath");
    auto sci_name = tags.find("sci_name");
    if (strcmp(node->name, "")) {
      out << "@SQ\tSN:" << node->name <<
          "\tLN:" << kFakeContigLength <<
          (gi == tags.cend() ? "" : "\tAS:" + (*gi).second) <<
          (fastapath == tags.cend() ? "" : "\tUR:" + (*fastapath).second) <<
          (sci_name == tags.cend() ? "" : "\tSP:" + (*sci_name).second) <<
          std::endl;
    }
  }
}

void ReadProcessor::print_sam_line(int32_t id, const std::string& suffix, std::ostream& out) {
  out << read_name_ << "\t";
  if (id == -1) {
    out << "4\t*\t0\t0\t*\t";
  } else {
    const knhx1_t* node = tree_.node_by_id(id);
    out << "0\t" << node->name << "\t1\t60\t" << coverage_cigar_ << "\t";
  }
  out << "*\t0\t0\t" << read_ << "\t" << qualities_ << "\t";
  if (id != -1) {
    out << "h1:i:" << best_hit_ << "\t";
    out << "c1:i:" << best_coverage_ << "\t";
    out << "hc:Z:" << hit_cigar_;
  }
  out << suffix << "\n";
}

void ReadProcessor::print_kraken_line(int32_t id, std::ostream& out) {
  if (simulate_lca_) {
    std::cerr << "simulate_lca_ is not implemented yet" << std::endl;
    exit(1);
  }
  if (id == -1) {
    out << "U\t";
  } else {
    out << "C\t";
  }
  out << read_name_ << "\t";
  if (id == -1) {
    out << "0\t";
  } else {
    const knhx1_t* node = tree_.node_by_id(id);
    out << node->name << "\t";
  }
  out << read_length_ << "\t";
  out << krakmers_ << "\n";
}

void ReadProcessor::clear() {
  clear_masks();
  matching_nodes_.clear();
  best_hit_nodes_.clear();
  best_coverage_nodes_.clear();
}

void ReadProcessor::parse_krakline(const std::string& krakline) {
  // std::cerr << "process krakline: " << krakline << std::endl;
  std::vector<std::string> parts = split(krakline, '\t');
  read_name_ = parts[1];
  read_length_ = std::stoi(parts[3]);
  krakmers_ = parts[4];
  if (parts.size() == 7) {
    read_ = parts[5];
    qualities_ = parts[6];
  } else {
    read_ = "*";
    qualities_ = "*";
  }
  fill_masks_from_kmer_blocks();
}

void ReadProcessor::propagate_matching_kmers() {
  for (auto id : matching_nodes_) {
    for (auto upper_node_id : tree_.upper_nodes(id)) {
      if (matching_nodes_.count(upper_node_id) > 0) {
        copy_masks(upper_node_id, id);
      }
    }
  }
}

void ReadProcessor::copy_masks(int32_t node_from, int32_t node_to) {
  // std::cerr << "copy node_from: " << node_from << " node_to: " << node_to << std::endl;
  auto& hit_mask_from = hit_masks_[node_from];
  auto& hit_mask_to = hit_masks_[node_to];
  for (size_t i = 0; i < hit_mask_size(); ++i) {
    hit_mask_to[i] |= hit_mask_from[i];
  }

  auto& coverage_mask_from = coverage_masks_[node_from];
  auto& coverage_mask_to = coverage_masks_[node_to];
  for (size_t i = 0; i < coverage_mask_size(); ++i) {
    coverage_mask_to[i] |= coverage_mask_from[i];
  }
}

void ReadProcessor::fill_masks_from_kmer_blocks() {
  size_t kmer_count = 0;
  size_t block_start = 0;
  size_t block_end = 0;
  while (block_start < krakmers_.size()) {
    while (block_end < krakmers_.size() && krakmers_[block_end] != ' ') {
      block_end++;
    }
    size_t kmer_in_block = fill_masks_from_block(
        std::experimental::string_view(krakmers_.c_str() + block_start,
        block_end - block_start), kmer_count);
    kmer_count += kmer_in_block;
    block_start = block_end + 1;
    block_end = block_start;
  }
  if (kmer_count + k_ - 1 != read_length_ && read_length_ >= k_) {
    std::cerr << "read length does not correspond to kmers blocks total length" << std::endl;
    exit(1);
  }
}

size_t ReadProcessor::fill_masks_from_block(
    const std::experimental::string_view& block, size_t processed_kmer_count) {
  size_t colon_position = 0;
  size_t node_start = 0;
  size_t node_end = 0;
  bool no_nodes = false;
  bool ambiguous_kmers = false;
  std::vector<int32_t> node_ids;
  while (node_start < block.size()) {
    while (node_end < block.size() && block[node_end] != ',') {
      if (block[node_end] == ':') {
        colon_position = node_end;
        break;
      }
      node_end++;
    }
    std::experimental::string_view node_name(
        block.data() + node_start, node_end - node_start);
    if (node_name == "0") {
      no_nodes = true;
      break;
    }
    if (node_name == "A") {
      ambiguous_kmers = true;
      break;
    }
    int32_t id = tree_.id_by_name(node_name);
    node_ids.push_back(id);
    matching_nodes_.insert(id);
    if (colon_position > 0) {
      break;
    }
    node_start = node_end + 1;
    node_end = node_start;
  }
  size_t count = std::stoi(block.data() + colon_position + 1);
  // auto block_split = split(block, ':');
  // size_t count = std::stoi(block_split[1]);
  // std::string& ids = block_split[0];
  // std::vector<int32_t> node_ids;
  // bool no_nodes = false;
  // bool ambiguous_kmers = false;
  // for (auto& node_name : split(ids, ',')) {
  //   if (node_name == "0") {
  //     no_nodes = true;
  //     break;
  //   }
  //   if (node_name == "A") {
  //     ambiguous_kmers = true;
  //     break;
  //   }
  //   int32_t id = tree_.id_by_name(node_name);
  //   node_ids.push_back(id);
  //   matching_nodes_.insert(id);
  // }
  if (!no_nodes && !ambiguous_kmers) {
    if (simulate_lca_) {
      // not implemented yet
      std::cerr << "simulate_lca is not implemented yet" << std::endl;
      exit(1);
    }
    set_masks(node_ids, processed_kmer_count, count);
  }
  return count;
  // print_masks();
}

void ReadProcessor::clear_masks() {
  for (auto id : matching_nodes_) {
    std::fill(hit_masks_[id].begin(), hit_masks_[id].end(), 0);
    std::fill(coverage_masks_[id].begin(), coverage_masks_[id].end(), 0);
  }
}

void ReadProcessor::print_masks() const {
  for (auto id : matching_nodes_) {
    auto& hit_mask = hit_masks_[id];
    auto& coverage_mask = coverage_masks_[id];
    std::cerr << tree_.node_by_id(id)->name << " " << hit_mask_size() << " " << coverage_mask_size() << std::endl;
    for (size_t i = 0; i < hit_mask_size(); ++i) {
      std::cerr << hit_mask[i];
    }
    std::cerr << std::endl;
    for (size_t i = 0; i < coverage_mask_size(); ++i) {
      std::cerr << coverage_mask[i];
    }
    std::cerr << std::endl;
  }
}

void ReadProcessor::set_masks(const std::vector<int32_t>& node_ids, size_t block_position, size_t block_length) {
  for (auto id : node_ids) {
    auto& hit_mask = hit_masks_[id];
    if (hit_mask.size() <= hit_mask_size()) {
      hit_mask.resize(2 * hit_mask_size(), 0);
    }
    for (size_t i = block_position; i < block_position + block_length; ++i) {
      hit_mask[i] = 1;
    }
    auto& coverage_mask = coverage_masks_[id];
    if (coverage_mask.size() <= coverage_mask_size()) {
      coverage_mask.resize(2 * coverage_mask_size(), 0);
    }
    for (size_t i = block_position; i < block_position + block_length + k_ - 1; ++i) {
      coverage_mask[i] = 1;
    }
  }
}

void ReadProcessor::fill_cigar(const std::vector<int16_t>& mask, size_t mask_size, std::string& cigar) {
  cigar.clear();
  size_t block_start = 0;
  while (block_start < mask_size) {
    size_t position = block_start;
    while (position < mask_size && mask[block_start] == mask[position]) {
      position++;
    }
    char match = mask[block_start] == 0 ? 'X' : '=';
    cigar += std::to_string(position - block_start);
    cigar += match;
    block_start = position;
  }
}
