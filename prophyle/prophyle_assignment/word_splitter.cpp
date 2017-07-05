#include "word_splitter.h"
#include <sstream>

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> parts;
  std::stringstream ss;
  ss.str(s);
  std::string part;
  while (std::getline(ss, part, delim)) {
    parts.push_back(part);
  }
  return parts;
}
