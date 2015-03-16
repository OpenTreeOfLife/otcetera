#if !defined OTCETERA_BASE_INCLUDES_H
#define OTCETERA_BASE_INCLUDES_H
#include <cassert>
#include <set>
#include <functional>
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored  "-Wweak-vtables"
#define ELPP_CUSTOM_COUT std::cerr
#define NOT_IMPLEMENTED assert("not implemented"[0] == 'f');
#include "otc/easylogging++.hpp"

namespace otc {
extern bool debuggingOutputEnabled;
extern long ottIDBeingDebugged;
// Might want to move to using
//	https://github.com/lczech/genesis/blob/master/src/tree/bipartitions.hpp
// at some point for faster operations on sets of indices

//using OttIdUSet = std::unordered_set<long>;
using OttId = long;
using OttIdOSet = std::set<long>;
using OttIdSet = OttIdOSet;

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
