#ifndef OTCETERA_WRITE_DOT_H
#define OTCETERA_WRITE_DOT_H

#include <map>
#include <vector>
#include <vector>
#include "otc/otc_base_includes.h"
namespace otc {

void writeDOTForEmbedding(std::ostream & out,
                          const NodeWithSplits * nd,
                          const std::vector<TreeMappedWithSplits *> &,
                          const std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> & eForNd,
                          bool entireSubtree,
                          bool includeLastTree);

} //namespace
#endif
