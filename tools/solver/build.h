#ifndef BUILD_H
#define BUILD_H

#include <vector>
#include "rsplit.h"
#include "solution.h"

bool BUILD(std::shared_ptr<Solution>& solution, const std::vector<ConstRSplit>& new_splits);
extern std::vector<int> indices;

#endif
