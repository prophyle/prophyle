#include "tree_index.h"
#include "read_processor.h"
#include <iostream>
#include <fstream>
#include <getopt.h>

struct Arguments {
  size_t k;
  std::string input_file = "-";
  std::string newick_file;
  AssignmentOutputFormat format = AssignmentOutputFormat::Kraken;
  Measure criteria = Measure::C1;
  bool simulate_lca = false;
  bool annotate = false;
  bool tie_lca = false;
  bool not_translate_blocks = false;
};

Arguments parse_arguments(int argc, char *argv[]) {
  int32_t c;
  int index;
  Arguments arguments;
  const struct option longopts[] =
      {
          {"input", required_argument, 0, 'i'},
          {"kmer-size", required_argument, 0, 'k'},
          {"newick-tree", required_argument, 0, 'n'},
          {"oformat", required_argument, 0, 'f'},
          {"measure", required_argument, 0, 'm'},
          {"sim-lca", optional_argument, 0, 'l'},
          {"annotate", optional_argument, 0, 'i'},
          {"tie-lca", optional_argument, 0, 't'},
          {"nontransl-blocks", optional_argument, 0, 'd'},
          {0,0,0,0},
      };
  while ((c = getopt_long(argc, argv, "i:k:n:f:m:latd", longopts, &index)) >= 0) {
    switch (c) {
      case 'i': arguments.input_file = optarg; break;
      case 'k': {
        int32_t k_int = std::stoi(optarg);
        if (k_int <= 0) {
          std::cerr << "k should be > 0" << std::endl;
          exit(1);
        }
        arguments.k = static_cast<size_t>(k_int);
        break;
      }
      case 'n': arguments.newick_file = optarg; break;
      case 'f': {
        std::string format_str(optarg);
        if (format_str == "sam") {
          arguments.format = AssignmentOutputFormat::Sam;
        } else if (format_str == "kraken") {
          arguments.format = AssignmentOutputFormat::Kraken;
        } else {
          std::cerr << "assignment format should be sam or kraken" << std::endl;
          exit(1);
        }
        break;
      }
      case 'm': {
        std::string measure_str(optarg);
        if (measure_str == "c1") {
          arguments.criteria = Measure::C1;
        } else if (measure_str == "h1") {
          arguments.criteria = Measure::H1;
        } else {
          std::cerr << "measure should be c1 or h1" << std::endl;
          exit(1);
        }
        break;
      }
      case 'l': arguments.simulate_lca = true; break;
      case 'a': arguments.annotate = true; break;
      case 't': arguments.tie_lca = true; break;
      case 'd': arguments.not_translate_blocks = true; break;
      default: {
        std::cerr << "argument " << c << " is not supported" << std::endl;
        exit(1);
      };
    }
  }
  return arguments;
}

int main(int argc, char *argv[]) {
  Arguments arguments = parse_arguments(argc, argv);
  TreeIndex tree = TreeIndex(arguments.newick_file);
  ReadProcessor read_processor = ReadProcessor(tree, arguments.k, arguments.simulate_lca, arguments.annotate,
      arguments.tie_lca, arguments.not_translate_blocks);

  if (arguments.input_file != "-") {
    freopen(arguments.input_file.c_str(), "r", stdin);
  }

  std::string line;
  while (getline(std::cin, line)) {
    read_processor.process_krakline(line, arguments.format, arguments.criteria);
  }

//  std::string krakline = "U\tread1\t0\t60\tn001:21 n001,n154:1 n001:8 n001,n410:1 n001:12 n001,n151:1 n001:1";
//  read_processor.process_krakline(krakline, arguments.format, arguments.criteria);
//  krakline = "U\tread2\t0\t37\tn001:3 n001,n023:2 n001:1 n001,n154,n176:1 n001,n176:2 n001,n154,n176:4 n001,n176:2 n001,n154,n176:4 n001,n176:1 n001:2";
//  read_processor.process_krakline(krakline, arguments.format, arguments.criteria);
//  krakline = "U\tread3\t0\t60\tn002:20 0:11 n006:7 n006,n134:1 n006:6";
//  read_processor.process_krakline(krakline, arguments.format, arguments.criteria);
//  krakline = "U\tread4\t0\t60\tn175,n199,n223,n271,n272,n287,n293,n307,n323,n329,n332,n348,n354,n355,n358:45";
//  read_processor.process_krakline(krakline, arguments.format, arguments.criteria);
//  krakline = "U\tread5\t0\t60\t0:45U";
//  read_processor.process_krakline(krakline, arguments.format, arguments.criteria);

  return 0;
}
