/*
  prophyle_index build command high-level implementation.
  Author: Kamil Salikhov <salikhov.kamil@gmail.com>
  Licence: MIT
*/

#ifndef PROPHYLE_INDEX_BUILD_H
#define PROPHYLE_INDEX_BUILD_H

#include "prophyle_utils.h"

void build_index(const char *prefix, const prophyle_index_opt_t *opt, int sa_intv);

#endif //PROPHYLE_INDEX_BUILD_H
