#if !defined OTCETERA_BASE_INCLUDES_H
#define OTCETERA_BASE_INCLUDES_H
#include "assert.hh"
#include <set>
#include <functional>
#ifdef __clang__
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored  "-Wweak-vtables"
#endif
#define ELPP_CUSTOM_COUT std::cerr
#define ELPP_THREAD_SAFE 1
#define ELPP_STACKTRACE_ON_CRASH 1

#include "otc/easylogging++.hpp"

#define OTC_UNREACHABLE {LOG(ERROR)<<"Unreachable code reached!"; std::abort();}

namespace otc {
extern bool debugging_output_enabled;


void throw_ott_id_type_too_small_exception(long);
#if defined(LONG_OTT_ID)
    using OttId = long;
    inline OttId check_ott_id_size(long raw_ott_id) {
        return raw_ott_id;
    }
#else
    using OttId = int;
    inline OttId check_ott_id_size(long raw_ott_id) {
        if (raw_ott_id >= std::numeric_limits<OttId>::max()) {
            throw_ott_id_type_too_small_exception(raw_ott_id);
        }
        return static_cast<OttId>(raw_ott_id);
    }
#endif
using OttIdSet = std::set<OttId>;

// forward decl
class RTSplits;
class RTNodeNoData;
class RTreeNoData;
template<typename T> class RootedTreeNode;
template<typename T> class RTreeOttIDMapping;
template<typename T, typename U> class RootedTree;
template<typename T> class RTreeOttIDMapping;
template<typename T, typename U> class NodePairing;
template<typename T, typename U> class PathPairing;
template<typename T, typename U> class NodeEmbedding;
template<typename T, typename U> class SupertreeContext;
template<typename T, typename U> class RootedForest;

struct RTNodeSmallestChild {
    OttId smallest_child = 0;
};

using RootedTreeNodeNoData = RootedTreeNode<RTNodeNoData>;
using RootedTreeTopologyNoData = RootedTree<RTNodeNoData, RTreeNoData> ;
using NodeWithSplits = RootedTreeNode<RTSplits>;
using NodeWithSplitsPred = std::function<bool(const NodeWithSplits &)>;
using MappedWithSplitsData = RTreeOttIDMapping<RTSplits>;
using MappedWithEmptyNodeData = RTreeOttIDMapping<RTNodeNoData>;
using TreeMappedWithSplits = RootedTree<RTSplits, MappedWithSplitsData>;
using TreeMappedEmptyNodes = RootedTree<RTNodeNoData, MappedWithEmptyNodeData> ;
using TreeMappedWithSplits = RootedTree<RTSplits, MappedWithSplitsData>;
using SupertreeContextWithSplits = SupertreeContext<NodeWithSplits, NodeWithSplits>;
using NodePairingWithSplits = NodePairing<NodeWithSplits, NodeWithSplits>;
using PathPairingWithSplits = PathPairing<NodeWithSplits, NodeWithSplits>;
using NodeEmbeddingWithSplits = NodeEmbedding<NodeWithSplits, NodeWithSplits>;

} // namespace otc
#endif
