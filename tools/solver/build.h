#ifndef BUILD_H
#define BUILD_H

#include <vector>
#include <optional>
#include "rsplit.h"
#include "solution.h"
#include "tree.h"

// Add new_splits to an existing solution
bool BUILDINC(std::shared_ptr<Solution>& solution, const std::vector<ConstRSplit>& new_splits);

// Run BUILD and report if the splits are consistent
bool BUILD_check(const std::vector<int> all_leaves_indices, const std::vector<ConstRSplit>& new_splits);

// Run BUILD and report if the splits are consistent
std::unique_ptr<Tree_t> BUILD(const std::vector<int> all_leaves_indices, const std::vector<ConstRSplit>& new_splits);

extern std::vector<int> indices;

#endif
