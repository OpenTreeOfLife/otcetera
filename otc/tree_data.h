#ifndef OTCETERA_TREE_DATA_H
#define OTCETERA_TREE_DATA_H
// Classes that can serve as the template args for trees and nodes
#include <map>
#include <set>
#include "otc/otc_base_includes.h"

namespace otc {
template<typename, typename> class RootedTree;
template<typename> class RootedTreeNode;

template<typename T>
class RTreeOttIDMapping {
    public:
        typedef RootedTreeNode<T> NodeType;
        std::map<OttId, NodeType *> ott_id_to_node;
        // if a node is pruned, the entry in ott_id_to_node will refer to an alias
        //   if the alias is later pruned, we need to know what nodes it is aliasing
        //   so that ott_id_to_node does not point to dangling nodes.
        std::map<NodeType *, std::set<OttId> > is_alias_for;
        std::map<OttId, NodeType *> ott_id_to_detached_node;
        bool des_id_sets_contain_internals;
        NodeType * get_node_by_ott_id(OttId ottId) const {
            const auto it = ott_id_to_node.find(ottId);
            return (it == ott_id_to_node.end() ? nullptr : it->second);
        }
};

class RTSplits {
    public:
    std::set<OttId> des_ids;
    int depth = 0;
};


template<typename T>
inline void verify_ott_id_mapping(const T & tree) {
    for (const auto & mp : tree.get_data().ott_id_to_node) {
        if (mp.first != mp.second->get_ott_id()) {
            assert(mp.first == mp.second->get_ott_id());
            throw OTCError("Cache of OTT ID->node is stale " + std::to_string(mp.first) + " != " + std::to_string(mp.second->get_ott_id()));
        }
    }
}

template<typename T>
inline std::vector<typename T::node_type *> get_nodes_aliased_by(typename T::node_type *nd, const T & tree) {
    const auto & iaf = tree.get_data().is_alias_for;
    const auto aIt = iaf.find(nd);
    if (aIt == iaf.end()) {
        std::vector<typename T::node_type *> empty;
        return empty;
    }
    std::vector<typename T::node_type *> r;
    r.reserve(aIt->second.size());
    const auto & o2d = tree.get_data().ott_id_to_detached_node;
    for (auto oi : aIt->second) {
        typename T::node_type * detached = o2d.find(oi)->second;
        r.push_back(detached);
    }
    return r;
}

} // namespace otc
#endif
