#ifndef SOLVER_TREE_H
#define SOLVER_TREE_H

#include "otc/supertree_util.h"
#include "otc/tree_iter.h"

typedef otc::TreeMappedWithSplits Tree_t;
typedef Tree_t::node_type node_t;

inline int depth(const Tree_t::node_type* nd)
{
    return nd->get_data().depth;
}

#endif
