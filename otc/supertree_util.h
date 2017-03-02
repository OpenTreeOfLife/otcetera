#ifndef OTCETERA_SUPERTREE_UTIL_H
#define OTCETERA_SUPERTREE_UTIL_H

#include <map>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/util.h"
#include "otc/tree_iter.h"
#include "otc/tree_data.h"
#include <boost/optional.hpp>

namespace otc {
std::unique_ptr<TreeMappedWithSplits> cloneTree(const TreeMappedWithSplits &);

template<typename T, typename U>
void updateAncestralPathOttIdSet(T * nd,
                                const OttIdSet & oldEls,
                                const OttIdSet & newEls,
                                std::map<const T *, NodeEmbedding<T, U> > & m);

template<typename T>
bool canBeResolvedToDisplayIncExcGroup(const T *nd, const OttIdSet & incGroup, const OttIdSet & excGroup);
template<typename T>
bool canBeResolvedToDisplayOnlyIncGroup(const T *nd, const OttIdSet & incGroup);

template<typename T, typename U>
void reportOnConflicting(std::ostream & out,
                         const std::string & prefix,
                         const T * scaffold,
                         const std::set<PathPairing<T, U> *> & exitPaths,
                         const OttIdSet & phyloLeafSet);

// takes 2 "includeGroups" from different PhyloStatements.
//  `culled` has been pruned down to the common leafSet, and the common leafSet is passed in as `leafSet
bool culledAndCompleteIncompatWRTLeafSet(const OttIdSet & culled, const OttIdSet & complete, const OttIdSet & leafSet);

template<typename T, typename U>
void copyStructureToResolvePolytomy(const T * srcPoly,
                                    U & destTree,
                                    typename U::node_type * destPoly,
                                    SupertreeContextWithSplits *);
// assumes that nd is the mrca of incGroup and excGroup IDs
template<typename T>
bool canBeResolvedToDisplay(const T *nd, const OttIdSet & incGroup, const OttIdSet & leafSet);


// assumes that nd is the mrca of incGroup and excGroup IDs
template<typename T>
inline bool canBeResolvedToDisplay(const T *nd, const OttIdSet & incGroup, const OttIdSet & leafSet) {
    const OttIdSet excGroup = set_difference_as_set(leafSet, incGroup);
    return canBeResolvedToDisplayIncExcGroup(nd, incGroup, excGroup);
}

template<typename T>
void addDesIdsToNdAndAnc(T * nd, const OttIdSet & oid) {
    nd->get_data().desIds.insert(begin(oid), end(oid));
    for (auto anc : iter_anc(*nd)) {
        anc->get_data().desIds.insert(begin(oid), end(oid));
    }
}

template<typename T>
void removeDesIdsToNdAndAnc(T * nd, const OttIdSet & oid) {
    nd->get_data().desIds = set_difference_as_set(nd->get_data().desIds, oid);
    for (auto anc : iter_anc(*nd)) {
        anc->get_data().desIds = set_difference_as_set(anc->get_data().desIds, oid);
    }
}

// used post-suppression of monotypic taxa to create a map
//  from the alias to the original OTT ID.
template<typename T>
inline std::map<long, long> generateIdRemapping(const T & tree) {
    const auto & id2ndMap = tree.get_data().ottIdToNode;
    std::map<long, long> r;
    for (const auto & idNdPair : id2ndMap) {
        const auto & inID = idNdPair.first;
        const auto outID = idNdPair.second->getOttId();
        if (outID != inID) {
            assert(!contains(r, inID));
            r[inID] = outID;
        }
    }
    return r;
}

//currently not copying names
inline std::unique_ptr<TreeMappedWithSplits> cloneTree(const TreeMappedWithSplits &tree) {
    TreeMappedWithSplits * rawTreePtr = new TreeMappedWithSplits();
    try {
        NodeWithSplits * newRoot = rawTreePtr->createRoot();
        auto r = tree.getRoot();
        assert(r->hasOttId());
        newRoot->setOttId(r->getOttId());
        std::map<const NodeWithSplits *, NodeWithSplits *> templateToNew;
        templateToNew[r]= newRoot;
        std::map<long, NodeWithSplits *> & newMap = rawTreePtr->get_data().ottIdToNode;
        rawTreePtr->get_data().desIdSetsContainInternals = tree.get_data().desIdSetsContainInternals;
        for (auto nd : iter_pre_const(tree)) {
            auto p = nd->getParent();
            if (p == nullptr) {
                continue;
            }
            auto t2nIt = templateToNew.find(p);
            assert(t2nIt != templateToNew.end());
            auto ntp = t2nIt->second;
            auto nn = rawTreePtr->createChild(ntp);
            assert(templateToNew.find(nd) == templateToNew.end());
            templateToNew[nd] = nn;
            if (nd->hasOttId()) {
                nn->setOttId(nd->getOttId());
                newMap[nd->getOttId()] = nn;
            } else {
                assert(false);
                throw OTCError("asserts false but not enabled");
            }
            nn->get_data().desIds = nd->get_data().desIds;
        }
    } catch (...) {
        delete rawTreePtr;
        throw;
    }
    return std::unique_ptr<TreeMappedWithSplits>(rawTreePtr);
}

template<typename T>
void sortChildOrderByLowestDesOttId(T *nd);

template<typename T>
void sortChildOrderByLowestDesOttId(T *deepest) {
    std::map<T *, long> node2Id;
    std::set<T *> internals;
    for (auto nd : iter_post_n(*deepest)) {
        if (nd->isTip()) {
            assert(nd->hasOttId());
            node2Id[nd] = nd->getOttId();
        } else {
            long lm = LONG_MAX;
            for (auto c : iter_child_const(*nd)) {
                auto coid = node2Id.at(const_cast<T *>(c));
                lm = std::min(lm, coid);
            }
            node2Id[nd] = lm;
            internals.insert(nd);
        }
    }
    for (auto nd : internals) {
        std::map<long, T *> id2child;
        for (auto c : iter_child(*nd)) {
            auto coid = node2Id.at(c);
            auto i2csize = id2child.size();
            id2child[coid] = c;
            assert(id2child.size() == 1 + i2csize); // assumes tip IDs are unique
        }
        assert(!id2child.empty());

        // Remove all the children - they are remembered in the map
        while(nd->hasChildren())
            nd->getFirstChild()->detachThisNode();

        // Add the children back in sorted order
        for(const auto& x: id2child)
            nd->addChild(x.second);
        assert(node2Id.at(nd->getFirstChild()) == node2Id.at(nd));
    }
}

// returns true if all of the children of nd which intersect with incGroup do NOT intersect w/ excGroup.
// NOTE: `nd` is assumed to be a common anc of all IDs in incGroup!
template<typename T>
inline bool canBeResolvedToDisplayIncExcGroup(const T *nd, const OttIdSet & incGroup, const OttIdSet & excGroup) {
    for (auto c : iter_child_const(*nd)) {
        if (haveIntersection(incGroup, c->get_data().desIds) && haveIntersection(excGroup, c->get_data().desIds)) {
            return false;
        }
    }
    return true;
}

// returns true if all of the children of nd which intersect with incGroup do NOT intersect w/ excGroup.
// NOTE: `nd` is assumed to be a common anc of all IDs in incGroup!
template<typename T>
inline bool canBeResolvedToDisplayOnlyIncGroup(const T *nd, const OttIdSet & incGroup) {
    for (auto c : iter_child_const(*nd)) {
        const auto & cdi = c->get_data().desIds;
        if (haveIntersection(incGroup, cdi) && (!isSubset(cdi, incGroup))) {
            return false;
        }
    }
    return true;
}

/// Functions below here are hackier in that they codify systems for embedding IDs in newick labels

template<typename T>
inline std::map<std::string, long> createIdsFromNames(const T & taxonomy);
template<typename T>
inline std::map<std::string, long> createIdsFromNamesFromTrees(const T& treeColl);
// awkward should merge the following 2 functions and use template specialization
template<typename T>
void setIdsFromNamesAndRefresh(T& tree, const std::map<std::string, long>& name_to_id);
template<typename T>
void setIdsFromNames(T& tree, const std::map<std::string, long>& name_to_id);
template<typename T>
void fillIdMapFromNames(const T & taxonomy, std::map<std::string, long> & name_to_id, long &nextId, bool allowRep);
std::string addOttId(const std::string & s, long id);
template<typename T>
inline void relabelNodesWithOttId(T& tree);
/// Create a mapping from name -> id.
/// throws exception for repeated names or unnamed tips
template<typename T>
inline std::map<std::string, long> createIdsFromNames(const T& taxonomy) {
    long id = 1;
    std::map<std::string, long> name_to_id;
    fillIdMapFromNames(taxonomy, name_to_id, id, false);
    return name_to_id;
}

/// Create a mapping from name -> id.
/// throws exception for unnamed tips - does NOT verify that a name only occurs once in a tree!
template<typename T>
inline std::map<std::string, long> createIdsFromNamesFromTrees(const T& treeColl) {
    long id = 1;
    std::map<std::string, long> name_to_id;
    for (const auto & tree : treeColl) {
        fillIdMapFromNames(*tree, name_to_id, id, true);
    }
    return name_to_id;
}


template<typename T>
inline void fillIdMapFromNames(const T & tree, std::map<std::string, long> & name_to_id, long & nextId, bool allowRep) {
    for(auto nd: iter_post_const(tree)) {
        if (nd->get_name().size()) {
            const std::string name = nd->get_name();
            auto it = name_to_id.find(name);
            if (it != name_to_id.end()) {
                if (not allowRep) {
                    throw OTCError()<<"Tip label '"<<name<<"' occurs twice!";
                }
            } else {
                name_to_id[name] = nextId++;
            }
        } else if (nd->isTip()) {
            throw OTCError()<<"tip has no label!";
        }
    }
}

/// Set ids on the tree based on the name
template<typename T>
inline void setIdsFromNamesAndRefresh(T& tree, const std::map<std::string,long> & name_to_id) {
    for(auto nd: iter_post(tree)){
        if (nd->get_name().size()) {
            const auto name = nd->get_name();
            const auto it = name_to_id.find(name);
            if (it == name_to_id.end()) {
                throw OTCError()<<"Can't find label '"<<name<<"' in taxonomy!";
            }
            const auto id = it->second;
            nd->setOttId(id);
            tree.get_data().ottIdToNode[id] = nd;
        } else if (nd->isTip()){
            throw OTCError()<<"Tree tip has no label!";
        }
    }
    clearAndfillDesIdSets(tree);
}

/// Set ids on the tree based on the name
template<typename T>
inline void setIdsFromNames(T& tree, const std::map<std::string,long> & name_to_id) {
    for(auto nd: iter_post(tree)){
        if (nd->get_name().size()) {
            const auto name = nd->get_name();
            const auto it = name_to_id.find(name);
            if (it == name_to_id.end()) {
                throw OTCError()<<"Can't find label '"<<name<<"' in taxonomy!";
            }
            const auto id = it->second;
            nd->setOttId(id);
        } else if (nd->isTip()){
            throw OTCError()<<"Tree tip has no label!";
        }
    }
}


inline std::string addOttId(const std::string & s, long id) {
    std::string tag = "ott" + std::to_string(id);
    if (not s.size()) {
        return tag;
    } else {
        return s + " " + tag;
    }
}

template<typename T>
inline void relabelNodesWithOttId(T& tree) {
    for(auto nd: iter_pre(tree)){
        if (nd->hasOttId()){
            nd->setName(addOttId(nd->get_name(),nd->getOttId()));
        }
    }
}

template<typename N>
std::size_t countChildren(const N * nd);
template<typename T>
std::size_t countLeaves(const T& tree);
template<typename T>
std::size_t countLeavesSubTree(const T* node);
template<typename N>
int mark(const N* node);
template<typename N>
int& mark(N* node);
template<typename N>
bool is_marked(const N* node, int bits);
template<typename N>
void set_mark(N* node, int bits);
template<typename N>
std::size_t countMarkedChildren(const N* nd, int bits);

template<typename T>
inline std::size_t countLeaves(const T& tree) {
    std::size_t count = 0U;
    for(auto nd: iter_post_const(tree)) {
        if (nd->isTip()) {
            count++;
        }
    }
    return count;
}

template<typename N>
inline std::size_t countLeavesSubTree(const N * tree) {
    std::size_t count = 0U;
    for(auto nd: iter_post_n_const(tree)) {
        if (nd->isTip()) {
            count++;
        }
    }
    return count;
}

template<typename N>
std::size_t countChildren(const N * nd) {
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
inline std::size_t countMarkedChildren(const N* nd, int bits) {
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
boost::optional<std::string> getSourceNodeName(const std::string& name);

} // namespace
#endif
