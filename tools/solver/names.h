#ifndef OTC_SOLVER_NAMES_H
#define OTC_SOLVER_NAMES_H

#include <vector>
#include <memory>                  // for unique_ptr
#include <set>
#include "solver/tree.h"

void add_root_and_tip_names(Tree_t& summary, Tree_t& taxonomy);
void add_names(Tree_t& summary, const std::vector<const Tree_t::node_type*>& compatible_taxa);

#endif
