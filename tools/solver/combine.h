#ifndef OTC_COMBINE_H
#define OTC_COMBINE_H

#include <vector>
#include <memory>                  // for unique_ptr
#include <set>
#include "solver/tree.h"
#include "otc/otc_base_includes.h" // for OttId
#include "otc/otcli.h"

std::unique_ptr<Tree_t> combine(std::vector<std::unique_ptr<Tree_t>>& trees, const std::set<otc::OttId>& incertae_sedis, boost::program_options::variables_map& args);

#endif
