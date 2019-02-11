#ifndef OTCETERA_SUPERTREE_UTIL_H
#define OTCETERA_SUPERTREE_UTIL_H

#include <map>
#include <set>
#include <limits>
#include "otc/otc_base_includes.h"
#include "otc/util.h"
#include "otc/tree_iter.h"
#include "otc/tree_data.h"
#include <optional>

namespace otc {
std::unique_ptr<TreeMappedWithSplits> clone_tree(const TreeMappedWithSplits &);

template<typename T, typename U>
void update_ancestral_path_ott_id_set(T * nd,
                                const OttIdSet & oldEls,
                                const OttIdSet & newEls,
                                std::map<const T *, NodeEmbedding<T, U> > & m);

template<typename T>
bool can_be_resolved_to_display_inc_exc_group(const T *nd, const OttIdSet & incGroup, const OttIdSet & excGroup);
template<typename T>
bool can_be_resolved_to_display_only_inc_exc_group(const T *nd, const OttIdSet & incGroup);

template<typename T, typename U>
void report_on_conflicting(std::ostream & out,
                         const std::string & prefix,
                         const T * scaffold,
                         const std::set<PathPairing<T, U> *> & exitPaths,
                         const OttIdSet & phyloLeafSet);

// takes 2 "includeGroups" from different PhyloStatements.
//  `culled` has been pruned down to the common leaf_set, and the common leaf_set is passed in as `leaf_set
bool culled_and_complete_incompat_wrt_leaf_set(const OttIdSet & culled, const OttIdSet & complete, const OttIdSet & leaf_set);

template<typename T, typename U>
void copy_structure_to_resolve_polytomy(const T * srcPoly,
                                    U & destTree,
                                    typename U::node_type * destPoly,
                                    SupertreeContextWithSplits *);
// assumes that nd is the mrca of incGroup and excGroup IDs
template<typename T>
bool can_be_resolved_to_display(const T *nd, const OttIdSet & incGroup, const OttIdSet & leaf_set);


// assumes that nd is the mrca of incGroup and excGroup IDs
template<typename T>
inline bool can_be_resolved_to_display(const T *nd, const OttIdSet & incGroup, const OttIdSet & leaf_set) {
    const OttIdSet excGroup = set_difference_as_set(leaf_set, incGroup);
    return can_be_resolved_to_display_inc_exc_group(nd, incGroup, excGroup);
}

template<typename T>
void add_des_ids_to_node_and_anc(T * nd, const OttIdSet & oid) {
    nd->get_data().des_ids.insert(begin(oid), end(oid));
    for (auto anc : iter_anc(*nd)) {
        anc->get_data().des_ids.insert(begin(oid), end(oid));
    }
}

template<typename T>
void remove_des_ids_from_node_and_anc(T * nd, const OttIdSet & oid) {
    nd->get_data().des_ids = set_difference_as_set(nd->get_data().des_ids, oid);
    for (auto anc : iter_anc(*nd)) {
        anc->get_data().des_ids = set_difference_as_set(anc->get_data().des_ids, oid);
    }
}

// used post-suppression of monotypic taxa to create a map
//  from the alias to the original OTT ID.
template<typename T>
inline std::map<OttId, OttId> generate_id_remapping(const T & tree) {
    const auto & id2ndMap = tree.get_data().ott_id_to_node;
    std::map<OttId, OttId> r;
    for (const auto & idNdPair : id2ndMap) {
        const auto & inID = idNdPair.first;
        const auto outID = idNdPair.second->get_ott_id();
        if (outID != inID) {
            assert(!contains(r, inID));
            r[inID] = outID;
        }
    }
    return r;
}

//currently not copying names
inline std::unique_ptr<TreeMappedWithSplits> clone_tree(const TreeMappedWithSplits &tree) {
    TreeMappedWithSplits * rawTreePtr = new TreeMappedWithSplits();
    try {
        NodeWithSplits * newRoot = rawTreePtr->create_root();
        auto r = tree.get_root();
        assert(r->has_ott_id());
        newRoot->set_ott_id(r->get_ott_id());
        std::map<const NodeWithSplits *, NodeWithSplits *> templateToNew;
        templateToNew[r]= newRoot;
        std::map<OttId, NodeWithSplits *> & newMap = rawTreePtr->get_data().ott_id_to_node;
        rawTreePtr->get_data().des_id_sets_contain_internals = tree.get_data().des_id_sets_contain_internals;
        for (auto nd : iter_pre_const(tree)) {
            auto p = nd->get_parent();
            if (p == nullptr) {
                continue;
            }
            auto t2nIt = templateToNew.find(p);
            assert(t2nIt != templateToNew.end());
            auto ntp = t2nIt->second;
            auto nn = rawTreePtr->create_child(ntp);
            assert(templateToNew.find(nd) == templateToNew.end());
            templateToNew[nd] = nn;
            if (nd->has_ott_id()) {
                nn->set_ott_id(nd->get_ott_id());
                newMap[nd->get_ott_id()] = nn;
            } else {
                assert(false);
                throw OTCError("asserts false but not enabled");
            }
            nn->get_data().des_ids = nd->get_data().des_ids;
        }
    } catch (...) {
        delete rawTreePtr;
        throw;
    }
    return std::unique_ptr<TreeMappedWithSplits>(rawTreePtr);
}

template<typename T>
void sort_children_by_lowest_des_ott_id(T *nd);

template<typename T>
void sort_children_by_lowest_des_ott_id(T *deepest) {
    std::map<T *, OttId> node2Id;
    std::set<T *> internals;
    for (auto nd : iter_post_n(*deepest)) {
        if (nd->is_tip()) {
            assert(nd->has_ott_id());
            node2Id[nd] = nd->get_ott_id();
        } else {
            OttId lm = std::numeric_limits<OttId>::max();
            for (auto c : iter_child_const(*nd)) {
                auto coid = node2Id.at(const_cast<T *>(c));
                lm = std::min(lm, coid);
            }
            node2Id[nd] = lm;
            internals.insert(nd);
        }
    }
    for (auto nd : internals) {
        std::map<OttId, T *> id2child;
        for (auto c : iter_child(*nd)) {
            auto coid = node2Id.at(c);
            auto i2csize = id2child.size();
            id2child[coid] = c;
            assert(id2child.size() == 1 + i2csize); // assumes tip IDs are unique
        }
        assert(!id2child.empty());
        // Remove all the children - they are remembered in the map
        while(nd->has_children()) {
            nd->get_first_child()->detach_this_node();
        }
        // Add the children back in sorted order
        for(const auto& x: id2child) {
            nd->add_child(x.second);
        }
        assert(node2Id.at(nd->get_first_child()) == node2Id.at(nd));
    }
}

// returns true if all of the children of nd which intersect with incGroup do NOT intersect w/ excGroup.
// NOTE: `nd` is assumed to be a common anc of all IDs in incGroup!
template<typename T>
inline bool can_be_resolved_to_display_inc_exc_group(const T *nd, const OttIdSet & incGroup, const OttIdSet & excGroup) {
    for (auto c : iter_child_const(*nd)) {
        if (have_intersection(incGroup, c->get_data().des_ids) && have_intersection(excGroup, c->get_data().des_ids)) {
            return false;
        }
    }
    return true;
}

// returns true if all of the children of nd which intersect with incGroup do NOT intersect w/ excGroup.
// NOTE: `nd` is assumed to be a common anc of all IDs in incGroup!
template<typename T>
inline bool can_be_resolved_to_display_only_inc_exc_group(const T *nd, const OttIdSet & incGroup) {
    for (auto c : iter_child_const(*nd)) {
        const auto & cdi = c->get_data().des_ids;
        if (have_intersection(incGroup, cdi) && (!is_subset(cdi, incGroup))) {
            return false;
        }
    }
    return true;
}

/// Functions below here are hackier in that they codify systems for embedding IDs in newick labels

template<typename T>
inline std::map<std::string, OttId> create_ids_from_names(const T & taxonomy);
template<typename T>
inline std::map<std::string, OttId> create_ids_from_names_from_trees(const T& treeColl);
// awkward should merge the following 2 functions and use template specialization
template<typename T>
void set_ids_from_names_and_refresh(T& tree, const std::map<std::string, OttId>& name_to_id);
template<typename T>
void set_ids_from_names(T& tree, const std::map<std::string, OttId>& name_to_id);
template<typename T>
void fill_id_map_from_names(const T & taxonomy, std::map<std::string, OttId> & name_to_id, OttId &nextId, bool allowRep);
std::string add_ott_id(const std::string & s, OttId id);
template<typename T>
inline void relabel_nodes_with_ott_id(T& tree);
/// Create a mapping from name -> id.
/// throws exception for repeated names or unnamed tips
template<typename T>
inline std::map<std::string, OttId> create_ids_from_names(const T& taxonomy) {
    OttId id = 1;
    std::map<std::string, OttId> name_to_id;
    fill_id_map_from_names(taxonomy, name_to_id, id, false);
    return name_to_id;
}

/// Create a mapping from name -> id.
/// throws exception for unnamed tips - does NOT verify that a name only occurs once in a tree!
template<typename T>
inline std::map<std::string, OttId> create_ids_from_names_from_trees(const T& treeColl) {
    OttId id = 1;
    std::map<std::string, OttId> name_to_id;
    for (const auto & tree : treeColl) {
        fill_id_map_from_names(*tree, name_to_id, id, true);
    }
    return name_to_id;
}


template<typename T>
inline void fill_id_map_from_names(const T & tree, std::map<std::string, OttId> & name_to_id, OttId & nextId, bool allowRep) {
    for(auto nd: iter_post_const(tree)) {
        if (nd->get_name().size()) {
            const std::string name = nd->get_name();
            auto it = name_to_id.find(name);
            if (it != name_to_id.end()) {
                if (not allowRep) {
                    throw OTCError() << "Tip label '" << name << "' occurs twice!";
                }
            } else {
                name_to_id[name] = nextId++;
            }
        } else if (nd->is_tip()) {
            throw OTCError() << "tip has no label!";
        }
    }
}

/// Set ids on the tree based on the name
template<typename T>
inline void set_ids_from_names_and_refresh(T& tree, const std::map<std::string,OttId> & name_to_id) {
    for(auto nd: iter_post(tree)){
        if (nd->get_name().size()) {
            const auto name = nd->get_name();
            const auto it = name_to_id.find(name);
            if (it == name_to_id.end()) {
                throw OTCError() << "Can't find label '" << name << "' in taxonomy!";
            }
            const auto id = it->second;
            nd->set_ott_id(id);
            tree.get_data().ott_id_to_node[id] = nd;
        } else if (nd->is_tip()){
            throw OTCError() << "Tree tip has no label!";
        }
    }
    clear_and_fill_des_ids(tree);
}

/// Set ids on the tree based on the name
template<typename T>
inline void set_ids_from_names(T& tree, const std::map<std::string,OttId> & name_to_id) {
    for(auto nd: iter_post(tree)){
        if (nd->get_name().size()) {
            const auto name = nd->get_name();
            const auto it = name_to_id.find(name);
            if (it == name_to_id.end()) {
                throw OTCError() << "Can't find label '" << name << "' in taxonomy!";
            }
            const auto id = it->second;
            nd->set_ott_id(id);
        } else if (nd->is_tip()){
            throw OTCError() << "Tree tip has no label!";
        }
    }
}


inline std::string add_ott_id(const std::string & s, OttId id) {
    std::string tag = "ott" + std::to_string(id);
    if (not s.size()) {
        return tag;
    } else {
        return s + " " + tag;
    }
}

template<typename T>
inline void relabel_nodes_with_ott_id(T& tree) {
    for(auto nd: iter_pre(tree)){
        if (nd->has_ott_id()){
            nd->set_name(add_ott_id(nd->get_name(),nd->get_ott_id()));
        }
    }
}

template<typename N>
std::size_t count_children(const N * nd);
template<typename T>
std::size_t count_leaves(const T& tree);
template<typename T>
std::size_t count_leaves_subtree(const T* node);
template<typename N>
int mark(const N* node);
template<typename N>
int& mark(N* node);
template<typename N>
bool is_marked(const N* node, int bits);
template<typename N>
void set_mark(N* node, int bits);
template<typename N>
std::size_t count_marked_children(const N* nd, int bits);

template<typename T>
inline std::size_t count_leaves(const T& tree) {
    std::size_t count = 0U;
    for(auto nd: iter_post_const(tree)) {
        if (nd->is_tip()) {
            count++;
        }
    }
    return count;
}

template<typename N>
inline std::size_t count_leaves_subtree(const N * tree) {
    std::size_t count = 0U;
    for(auto nd: iter_post_n_const(tree)) {
        if (nd->is_tip()) {
            count++;
        }
    }
    return count;
}

template<typename N>
std::size_t count_children(const N * nd) {
    std::size_t count = 0;
    for(auto nd2: iter_child_const(*nd)){
        count++;
    }
    return count;
}

template<typename N>
inline int mark(const N* node) {
    return node->get_data().mark;
}

template<typename N>
inline int& mark(N* node) {
    return node->get_data().mark;
}

template<typename N>
inline bool is_marked(const N* node, int bits) {
    return (mark(node)&bits) == bits;
}

template<typename N>
inline void set_mark(N* node, int bits) {
    mark(node) |= bits;
}

template<typename N>
inline std::size_t count_marked_children(const N* nd, int bits) {
    std::size_t count = 0;
    for(auto nd2: iter_child_const(*nd)) {
        if (is_marked(nd2, bits)) {
            count++;
        }
    }
    return count;
}

std::string study_from_tree_name(const std::string& name);
std::string tree_in_study_from_tree_name(const std::string& name);
std::string string_between_chars(const std::string & s, char beforeC, char endC);
std::string source_from_tree_name(const std::string & name);
std::optional<std::string> get_source_node_name(const std::string& name);

} // namespace
#endif
