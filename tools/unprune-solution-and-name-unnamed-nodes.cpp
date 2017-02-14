#include <algorithm>
#include <unordered_set>
#include <iterator>
#include <sstream>
#include <stack>
#include <tuple>
#include <deque>
#include "json.hpp"
#include "otc/otcli.h"
#include "otc/node_naming.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"

using namespace otc;

static std::string lostTaxaJSONFilename;
static std::string statsJSONFilename;
static std::string incertaeSedisFilename;

struct RTNodePartialDesSet {
    std::set<long> desIds;
    OttIdSet * nonexcluded_ids = nullptr; // union of IDs that descend from any child of an ancestor marked as incertae sedis
    long smallestChild  = 0;
};

struct RTTreeIncertaeSedisHolder {
    std::list<OttIdSet> ownedIdSets; // used for memory management.
};

using Node_t = RootedTreeNode<RTNodePartialDesSet>;
using Tree_t = RootedTree<RTNodePartialDesSet, RTTreeIncertaeSedisHolder>;

inline long smallestChild(const Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

inline long& smallestChild(Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

// Stores pointers for a taxon X which not in the solution.
//   points to the MRCA, and all of the attachment points for 
//   subtrees that are part of this taxon.
// An attachment point is:
//      1. an ancestor of at least one member of the taxon X,
//      2. an ancestor of at least one taxon that is not contained in the taxon X.
//      3. has at least one child (an attached node) that consists entirely
//          of member of taxon X
class LostTaxonLocation {
    public:
    LostTaxonLocation(long oid=0, Tree_t::node_type *mrcaNd=nullptr)
        :ottId(oid),
        mrca(mrcaNd) {
    }
    void addAttachmentPoint(Tree_t::node_type *attach, Tree_t::node_type *child) {
        attachNode2AttachedVec[attach].push_back(child);
    }
    public: //@TMP should be private.
    long ottId = 0L;
    Tree_t::node_type * mrca = nullptr;
    using node_vec = std::vector<Tree_t::node_type *>;
    using attach2node_vec = std::map<Tree_t::node_type *, node_vec>;
    // for each attachment point, we store which subtrees for this taxon attach there.
    attach2node_vec attachNode2AttachedVec;

};
typedef std::map<int, LostTaxonLocation> LostTaxonMap;

template<typename T>
LostTaxonMap unpruneTaxa(T & taxonomy,
                         T & solution,
                         std::ostream * statsStreamPtr,
                         const OttIdSet & incertae_sedis_ids);

void writeLostTaxa(std::ostream & out, const LostTaxonMap & ltm);

// impl.

using json = nlohmann::json;
using std::deque;
using std::list;
using std::map;
using std::ostream;
using std::pair;
using std::set;
using std::size_t;
using std::stack;
using std::string;
using std::tuple;
using std::unique_ptr;
using std::vector;


void writeLostTaxa(ostream & out, const LostTaxonMap & ltm) {
    json lostTaxa;
    for (auto ltmEl : ltm) {
        const auto & ottId = ltmEl.first;
        const auto & ltl = ltmEl.second; 
        json lostTaxonLocation;
        lostTaxonLocation["mrca"] = ltl.mrca->getName();
        const auto attachmentPoints = ltl.attachNode2AttachedVec;
        json apjson;
        for (auto attachPoint: attachmentPoints) {
            const auto & parent = attachPoint.first;
            const auto & childVec = attachPoint.second;
            json childVecJSON = json::array();
            for (auto c : childVec) {
                childVecJSON.push_back(c->getName());
            }
            apjson[parent->getName()] = childVecJSON;
        }
        lostTaxonLocation["attachment_points"] = apjson;
        lostTaxa["ott" + std::to_string(ottId)] = lostTaxonLocation;
    }
    json document;
    document["non_monophyletic_taxa"] = lostTaxa;
    out << document.dump(1) << std::endl;
}

// adds ottId to the desIds field of every node from firstNd down to ancAndLast (inclusive)
template <typename N>
void addToDesIdsForAnc(long ottId, N *firstNd, N *ancAndLast) {
    assert(firstNd);
    while (true) {
        firstNd->getData().desIds.insert(ottId);
        //LOG(DEBUG) << "adding " << ottId << " to node " << firstNd->getName(); dbWriteOttSet(" its desIds", firstNd->getData().desIds);
        if (firstNd == ancAndLast) {
            return;
        }

	// If this triggers, then ancAndLast is not an ancestor of firstNd
	assert(firstNd->getParent());

        firstNd = firstNd->getParent();
    }
}

template <typename N>
vector<N *> higherTaxPreOrderBelowBoundaries(N * root,
                                                  const set<N *> & boundaries) {
    set<N *> seen;
    vector<N *> r;
    N * curr = root;
    assert(!curr->getData().desIds.empty());
    stack<N *> toDealWith;
    while (true) {
        if (seen.count(curr) > 0) {
            if (toDealWith.empty()) {
                return r;
            }
            curr = toDealWith.top();
            toDealWith.pop();
            assert (seen.count(curr) == 0);
        } else {
            seen.insert(curr);
            if (boundaries.count(curr) == 0) {
                //LOG(DEBUG) << "higherTaxPreOrderBelowBoundaries adding " << curr->getOttId();
                r.push_back(curr);
                for (auto c : iter_child(*curr)) {
                    if (!c->getData().desIds.empty()) {
                        toDealWith.push(c);
                    }
                }
            }
        }
    }
}

template <typename N>
N * introduceMonotypicParent(N * taxon, N * solnChild) {
    assert(taxon);
    assert(solnChild);
    N * nn = new N(nullptr);
    solnChild->replaceThisNode(nn);
    nn->setName(taxon->getName());
    nn->setOttId(taxon->getOttId());
    nn->addChild(solnChild);
    nn->getData().desIds = solnChild->getData().desIds;
    return nn;
}

template <typename N>
N * addParentAndMoveUnsampledTaxChildren(N * taxon, N * solnChild) {
    auto * nn = introduceMonotypicParent(taxon, solnChild);
    moveUnsampledChildren(taxon, nn);
    return nn;
}

// Moves all of the children of taxon that have empty desIds fields
//  to be children of solnNode.
template <typename N>
void moveUnsampledChildren(N * taxon, N * solnNode) {
    assert(taxon);
    assert(solnNode);
    const auto children = all_children(taxon);
    for (auto child : children) {
        if (child->getData().desIds.empty()) {
            child->detachThisNode();
            solnNode->addChild(child);
        }
    }
}

// Walk tipward from the @mrca until we find subtrees whose
// leaves are purely descendants of @taxon.  Record these
// subtrees as the attachment points in @ltl.
template <typename N>
void registerAttachmentPoints(const N * taxon,
                              N * mrca,
                              LostTaxonLocation & ltl) {
    assert(taxon);
    assert(mrca);
    set<N *> visited;
    stack<N *> toDealWith;
    const auto & taxDes = taxon->getData().desIds;
    assert(taxDes.size() > 1);
    assert(isProperSubset(taxDes, mrca->getData().desIds));
    N * currSolnNd = mrca;
    // the root of the solution, must have all of the constituent taxa
    while (true) {
        if (visited.count(currSolnNd) > 0) {
            if (toDealWith.empty()) {
                return ;
            }
            currSolnNd = toDealWith.top();
            toDealWith.pop();
            assert (visited.count(currSolnNd) == 0);
        } else {
            visited.insert(currSolnNd);
            const auto & solDes = currSolnNd->getData().desIds;
            assert(solDes != taxDes); // we should not call this if a des of mrca = taxon
            for (auto c : iter_child(*currSolnNd)) {
                const auto & cdesId = c->getData().desIds;
                if (!areDisjoint(cdesId, taxDes)) {
                    if (isProperSubset(cdesId, taxDes)) {
                        ltl.addAttachmentPoint(currSolnNd, c);
                    } else {
                        toDealWith.push(c);
                    }
                }
            }
        }
    }
}

template <typename N>
void addChildrenOfNonMonophyleticTaxon(N * taxon,
                                       N * mrca,
                                       LostTaxonMap & ltm) {
    assert(taxon);
    assert(mrca);
    assert(taxon->hasOttId());
    const auto ottId = taxon->getOttId();
    assert(ltm.count(ottId) == 0); // should be adding a new element
    ltm[ottId] = LostTaxonLocation(ottId, mrca);
    auto & ltl = ltm[ottId];
    registerAttachmentPoints(taxon, mrca, ltl);
    moveUnsampledChildren(taxon, mrca);
}


template <typename N>
N * special_bisect_with_new_child(N * taxNode, N * currSolnNd, bool moveUnsamp) {
    N * nt = bisect_branch_with_new_child(currSolnNd);
    nt->setName(taxNode->getName());
    nt->setOttId(taxNode->getOttId());
    nt->getData().desIds = currSolnNd->getData().desIds;
    if (moveUnsamp) {
        moveUnsampledChildren(taxNode, nt);
    }
    return nt;
}

// returns the number of nodes on the backbone that are merged with this higher taxon
//      because the taxon was compatible (the same as) the node - will be either 0 or 1.
template <typename N>
size_t incorporateHigherTaxonNode(N* higherTaxonNd,
                                N* rootSolnNd,
                                set<N *> & nodesAddedForTaxa,
                                LostTaxonMap & ltm) {
    assert(higherTaxonNd);
    assert(rootSolnNd);
    LOG(DEBUG) << "higherTaxonNd->name = " << higherTaxonNd->getName();
    const auto & taxDes = higherTaxonNd->getData().desIds;
    assert(taxDes.size() > 0);
    N * currSolnNd = rootSolnNd;
    N * tipmostMRCA = nullptr;
    // the root of the solution, must have all of the constituent taxa
    assert(isSubset(taxDes, rootSolnNd->getData().desIds));
    while (true) {
        const auto & solDes = currSolnNd->getData().desIds;
        if (solDes == taxDes) {
            N * nt = nullptr;
            if (higherTaxonNd->isOutDegreeOneNode()) {
                if (currSolnNd == rootSolnNd) {
                    LOG(DEBUG) << "adding Monotypic child";
                    nt = special_bisect_with_new_child(higherTaxonNd, currSolnNd, false);
                } else {
                    LOG(DEBUG) << "introduceMonotypicParent";
                    nt = introduceMonotypicParent(higherTaxonNd, currSolnNd);
                }
            } else if (currSolnNd->hasOttId() 
                      && currSolnNd->getOttId() != higherTaxonNd->getOttId()) {
                if (currSolnNd == rootSolnNd) {
                    LOG(DEBUG) << "Adding forking child to make parent monotypic";
                    nt = special_bisect_with_new_child(higherTaxonNd, currSolnNd, true);
                } else {
                    LOG(DEBUG) << "addParentAndMoveUnsampledTaxChildren";
                    nt = addParentAndMoveUnsampledTaxChildren(higherTaxonNd, currSolnNd);
                }
            } else {
                LOG(DEBUG) << "moveUnsampledChildren";
                moveUnsampledChildren(higherTaxonNd, currSolnNd);
            }
            if (nt != nullptr) {
                nodesAddedForTaxa.insert(nt);
                return 0U;
            }
            return 1U;
        }
        assert(isProperSubset(taxDes, solDes));
        // look in this subtree for the higher taxon 
        tipmostMRCA = currSolnNd;
        N * nextSolnNd = nullptr;
        for (auto c : iter_child(*currSolnNd)) {
            const auto & cdesId = c->getData().desIds;
            //dbWriteOttSet("  c desIds", cdesId);        
            //dbWriteOttSet("  taxDes  ", taxDes);        
            if (!areDisjoint(cdesId, taxDes)) {
                //LOG(DEBUG) << "not disjoint nextSolnNd=" << (long) nextSolnNd;
                if (nextSolnNd == nullptr) {
                    nextSolnNd = c;
                } else {
                    // More than one child of currSolnNd, has a member of the taxon,
                    // so this is the mrca
                    // The fact that we are here (not in the branch above where we
                    //  tested for equality), means that this solution will not
                    //  display this taxon
                    assert(tipmostMRCA != nullptr);
                    addChildrenOfNonMonophyleticTaxon(higherTaxonNd, tipmostMRCA, ltm);
                    return 0U;
                }
            }
        }
        assert(nextSolnNd != nullptr);
        currSolnNd = nextSolnNd;
    }
}

using NumAddedTuple = tuple<size_t, size_t, size_t>;
// Moves unsampled taxa from the taxanomy (by finding them in `ott2tax`) to a slice of the
//  solution tree that is rooted at `rootSolnNd`.
// Uses the `rootSolnNd->getData().desIds` to denote the set of leaves of this slice of the
//  tree.
// If the taxon conflicts with the solution, then the unsampled taxa are attached at the
//  MRCA of the taxon in the solution, and an entry is added mapping the taxon's OTT Id
//  to a LostTaxonLocation object that explains where the elements of the taxon are located.
// returns a triple of numbers:
//      0 is number of new nodes added along the soln backbone (this does not include)
//          the internal nodes transferred from the taxonomy to the solution.
//      1 is the number of internal nodes that were merged with an existing node on the
//          backbone
//      2 is the number of solution leaves that were expanded to internals
template <typename N>
NumAddedTuple unpruneTaxaForSubtree(N *rootSolnNd,
                           const map<long, N*> & ott2tax, 
                           map<long, N*> & ott2soln, 
                           LostTaxonMap & ltm) {
    assert(rootSolnNd);
    assert(rootSolnNd->hasOttId());
    assert(!rootSolnNd->isTip());
    const auto ottId = rootSolnNd->getOttId();
    LOG(DEBUG) << "unpruneTaxaForSubtree for " << rootSolnNd->getName();
    assert(ott2soln.at(ottId) == rootSolnNd);
    N * rootTaxonNd = ott2tax.at(ottId);
    //LOG(DEBUG) << " root taxon is " << rootTaxonNd->getName();
    assert(!rootTaxonNd->isTip());
    // this will have the IDs for all of the include taxa for this slice of the tree
    auto & solnDesIds = rootSolnNd->getData().desIds;
    assert(!solnDesIds.empty());
    //dbWriteOttSet("solnRoot desIds    = ", solnDesIds);
    // desIds fields of the solution tree are filled in by the caller before this function.
    //  here we fill in those fields for the taxonomy nodes.
    //Here we add the sampled IDs to desId fields ofthe relevant taxa. This is potentially
    //  confusing because getData().desIds will hold only the IDs of taxa that are present in the soln
    //  tree (even when we are talking about the taxonomy node's desIds)
    // While we are walking through the "leaf" Ids for this tree slice, we'll also
    //  collect the leaf sets 
    set<N *> solnLeaves;
    set<N *> taxaLeaves;
    for (auto effectiveTipOttId : solnDesIds) {
        auto effTipTaxonNd = ott2tax.at(effectiveTipOttId);
        effTipTaxonNd->getData().desIds.clear();
        addToDesIdsForAnc(effectiveTipOttId,
                          effTipTaxonNd,
                          rootTaxonNd);
        solnLeaves.insert(ott2soln.at(effectiveTipOttId));
        taxaLeaves.insert(ott2tax.at(effectiveTipOttId));
    }
    //dbWriteOttSet("rootTaxonNd desIds = ", rootTaxonNd->getData().desIds);
    size_t numExpanded = 0;
    // If a tip of the solution is a higher taxon, then we should
    //  graft on the other tips here...
    for (auto l : solnLeaves) {
        if (l->isTip()) {
            assert(l->hasOttId());
            auto leafOttId = l->getOttId();
            auto taxonForLeaf = ott2tax.at(leafOttId);
            if (!taxonForLeaf->isTip()) {
                numExpanded += 1;
                auto cvec = all_children(taxonForLeaf);
                for (auto c : cvec) {
                    taxonForLeaf->removeChild(c);
                    l->addChild(c);
                }
            }
        }
    }
    size_t numMerged = numExpanded;
    // Walk through the taxonomy in a preorder fashion, but only accumulate
    //  a vector of those nodes with at least 1 member of the `desIds` field.
    //  these are the internal nodes that are induced by the taxa included
    //  in the solution tree.
    // We will reverse the vector in the next line, so that we can walk through
    //  in postorder.
    auto postOrderInTaxNd = higherTaxPreOrderBelowBoundaries(rootTaxonNd,
                                                             taxaLeaves);
    std::reverse(postOrderInTaxNd.begin(), postOrderInTaxNd.end());
    // Here we will add the relevant higher (non-leaf) taxa from the taxonomy to the 
    //  solution. It is crucial that we do this in postorder because if we have a series
    //  of nodes to introduce along a branch, the incorporateHigherTaxonNode function
    //  will add them below the attachment node (so postorder will assure that they 
    //  are correctly added in the tip->root orientation).
    set<N *> nodesAddedForTaxa;
    for (auto higherTaxonNd : postOrderInTaxNd) {
        if (higherTaxonNd == rootTaxonNd) {
            moveUnsampledChildren(higherTaxonNd, rootSolnNd);
            numMerged += 1 ;
        } else {
            numMerged += incorporateHigherTaxonNode(higherTaxonNd, rootSolnNd, nodesAddedForTaxa, ltm);
        }
    }
    // for the sake of a low memory footprint, here we 
    //  clear out the desIds fields for this slice of the tree.
    for (auto n : taxaLeaves) {
        n->getData().desIds.clear();
    }
    for (auto n : postOrderInTaxNd) {
        n->getData().desIds.clear();
    }
    set<N *> visited;
    for (auto n : solnLeaves) {
        auto c = n;
        while (visited.count(c) == 0 && c != rootSolnNd) {
            c->getData().desIds.clear();
            if (c != n) {
                visited.insert(c);
            }
            c = c->getParent();
        }
    }
    rootSolnNd->getData().desIds.clear();
    rootSolnNd->getData().desIds.insert(ottId);
    return NumAddedTuple(nodesAddedForTaxa.size(), numMerged, numExpanded);
}


struct UnpruneStats {
    size_t startNumNamedInternalsInSoln = 0U;
    size_t numUnnamedInternalsInSoln = 0U;
    size_t startNumSolnLeaves = 0U;
    size_t startNumSolnMonotypicInternals = 0U;
    size_t startNumSolnForkingInternals = 0U;
    size_t numInternalsAddedToBackbone = 0U;
    size_t numInternalsMergedOnBackbone = 0U;
    size_t numLeavesExpanded = 0U;
    size_t out_degree_many1 = 0U;
    size_t numTaxaLeaves = 0U;
    size_t numTaxaForkingInternals = 0U;
    size_t numTipsDesFromIncertaeSedis = 0U;
    size_t numInternalsDesFromIncertaeSedis = 0U;
    // maps inc_sed taxa nodes to tips in solution that are descendants of that taxon.
    map<Node_t *, set<Node_t *> > inc_sed_taxon_to_sampled_tips;
    // maps inc_sed taxon that is not monophyletic to MRCA in solution.
    map<Node_t *, Node_t *> non_monophyletic_inc_sed;
};

void fillNonexcludeIDFields(Tree_t & taxonomy, const OttIdSet & incertae_sedis_ids) {
    // Assuming that incertae sedis taxa are fairly rare, there will be lots of nodes that 
    //    share the same set of non-excluded taxa. Thus we have an owning store and lots
    //    of pointers to them.
    auto & owningList = taxonomy.getData().ownedIdSets;
    auto taxonomy_root_ptr = taxonomy.getRoot();
    owningList.push_back(OttIdSet());
    taxonomy_root_ptr->getData().nonexcluded_ids = &(*owningList.rbegin());
    for (auto nd: iter_pre(taxonomy)) {
        if (nd->isTip()) {
            continue;
        }
        auto & ndd = nd->getData();
        auto * curr_nonexc = ndd.nonexcluded_ids;
        assert(curr_nonexc != nullptr);
        set<Node_t *> inc_sed_children;
        set<Node_t *> normal_children;
        for (auto c : iter_child(*nd)) {
            if (incertae_sedis_ids.count(c->getOttId())) {
                inc_sed_children.insert(c);
            } else {
                normal_children.insert(c);
            }
        }
        if (inc_sed_children.empty()) {
            // all children can share same nonexcluded set as curr nd;
            for (auto c: normal_children) {
                c->getData().nonexcluded_ids = curr_nonexc;
            }
        } else {
            if (inc_sed_children.size() == 1) {
                // if there is only one, the incertae sedis child
                //     actually just has the same non-excluded as parent
                auto inc_sed_child = *inc_sed_children.begin();
                auto & inc_sed_child_data = inc_sed_child->getData();
                inc_sed_child_data.nonexcluded_ids = curr_nonexc;
            } else {
                // each incertae sedis child will have its own non-excluded set
                for (auto isc: inc_sed_children) {
                    owningList.push_back(*curr_nonexc);
                    isc->getData().nonexcluded_ids = &(*owningList.rbegin());
                }
                for (auto isc: inc_sed_children) {
                    const auto & isc_di = isc->getData().desIds;
                    for (auto oisc: inc_sed_children) {
                        if (oisc == isc) {
                            continue;
                        }
                        oisc->getData().nonexcluded_ids->insert(isc_di.begin(), isc_di.end());
                    }
                }
            }
            // all "normal" children share non-excluded, but it is larger than current node's
            if (!normal_children.empty()) {
                owningList.push_back(*curr_nonexc);
                auto normal_child_non_exc = &(*owningList.rbegin());
                for (auto isc: inc_sed_children) {
                    auto isc_di = isc->getData().desIds;
                    normal_child_non_exc->insert(isc_di.begin(), isc_di.end());
                }
                for (auto c: normal_children) {
                    c->getData().nonexcluded_ids = normal_child_non_exc;
                }
            }
        }
    }
}

void fillDesIdsForIncertaeSedisOnly(const map<long, Node_t *> & ott_to_tax,
                                    const OttIdSet & incertae_sedis_ids, 
                                    UnpruneStats & unprune_stats) {
    // As a part of implementing the incertae sedis logic, we also fill the desIds of the taxonomy.
    // This info is used to create the "non-excluded" sets, but it is expensive to do on the whole 
    //    taxonomy, so we just do it for subtrees rooted at an incertae sedis taxon...
    unprune_stats.numTipsDesFromIncertaeSedis = 0;
    unprune_stats.numInternalsDesFromIncertaeSedis = 0;
    for (auto ist_id: incertae_sedis_ids) {
        if (ott_to_tax.count(ist_id) == 0) {
            throw OTCError() << "Incertae sedis ID (" << ist_id << ") not in taxonomy.";
        }
        auto taxon = ott_to_tax.at(ist_id);
        for (auto nd: iter_post_n(*taxon)) {
            auto & di_ref = nd->getData().desIds;
            // empty check is just an optimization to avoid recalculating the results for a node.
            if (di_ref.empty()) {
                di_ref.insert(nd->getOttId());
                if (nd->isTip()) {
                    unprune_stats.numTipsDesFromIncertaeSedis += 1;
                } else {
                    unprune_stats.numInternalsDesFromIncertaeSedis += 1;
                    for (auto c : iter_child(*nd)) {
                        const auto & cdi_ref = c->getData().desIds;
                        di_ref.insert(cdi_ref.begin(), cdi_ref.end());
                    }
                }
            }
        }
    }
}


void indexNodesByOttId(Tree_t & taxonomy,
                       map<long, Node_t *> & ott_to_tax,
                       OttIdSet & monotypicOttIds,
                       UnpruneStats & unprune_stats) {
    for (auto nd: iter_post(taxonomy)){
        if (nd->hasOttId()) {
            const OttId taxon_ott_id = nd->getOttId();
            ott_to_tax[taxon_ott_id] = nd;
            if (nd->isOutDegreeOneNode()) {
                monotypicOttIds.insert(nd->getOttId());
            } else if (nd->isTip()) {
                unprune_stats.numTaxaLeaves += 1;
            } else {
                unprune_stats.numTaxaForkingInternals += 1;
            }
        } else {
            throw OTCError() << "There is a node in taxonomy without an OTT ID.\n";
        }
    }
}

// uses non-empty desId in tax nodes to figure out which descendants of any incertae sedis
//    taxon are sampled in the solution
void findIncertaeSedisInSolution(Tree_t & taxonomy,
                                 Tree_t & solution,
                                 const map<long, Node_t *> & ott_to_tax,
                                 UnpruneStats & unprune_stats) {
    auto & to_tips_map = unprune_stats.inc_sed_taxon_to_sampled_tips;
    for (auto snd: iter_leaf(solution)) {
        if (!snd->hasOttId()) {
            throw OTCError() << "Tip "<< snd->getName() << " in solution lacks an OTT ID";
        }
        auto ott_id = snd->getOttId();
        auto taxon_node = ott_to_tax.at(ott_id);
        if (taxon_node->getData().desIds.empty()) {
            continue;    
        }
        to_tips_map[taxon_node].insert(snd);
        auto anc = taxon_node->getParent();
        while (anc && !(anc->getData().desIds.empty())) {
            to_tips_map[anc].insert(snd);
            anc = anc->getParent();
        }
    }
}

template <typename N>
N * find_mrca_without_desids(const set<N *> & tip_node_set) {
    if (tip_node_set.size() < 2)  {
        if (tip_node_set.empty()) {
            return nullptr;
        }
        return *tip_node_set.begin();
    }
    auto nit = tip_node_set.begin();
    N * curr = *nit++;
    deque<N *> root_to_mrca;
    set<N *> seen_nodes;
    set<N *> dequed_nodes;
    while (curr != nullptr) {
        root_to_mrca.push_front(curr);
        seen_nodes.insert(curr);
        dequed_nodes.insert(curr);
        curr = curr->getParent();
    }
    for (; nit != tip_node_set.end(); ++nit) {
        curr = *nit;
        while (curr != nullptr) {
            if (seen_nodes.count(curr) != 0) {
                if (dequed_nodes.count(curr) != 0) {
                    while (root_to_mrca.back() != curr) {
                        auto td = root_to_mrca.back();
                        root_to_mrca.pop_back();
                        dequed_nodes.erase(td);
                        if (root_to_mrca.size() == 1) {
                            return root_to_mrca.back();
                        }
                    }
                }
                break;
            } else {
                seen_nodes.insert(curr);
            }
        }

    }
    return root_to_mrca.back();
}

// Taxa that are descendants of an incertae sedis taxon, will demand some special handling..
// so we find out which ones are not found in the solution before calling unprunePreppedInputs
// results are stored in unprune_stats. This relies on having been filled...
void findBrokenIncertaeSedisDescendants(Tree_t & taxonomy,
                                        Tree_t & solution,
                                        const map<long, Node_t *> & ott_to_tax,
                                        UnpruneStats & unprune_stats) {
    auto & to_tips_map = unprune_stats.inc_sed_taxon_to_sampled_tips;
    for (auto ttm_it: to_tips_map) {
        auto inc_sed_taxon = ttm_it.first;
        const auto & sampled = ttm_it.second;
        if (sampled.size() < 2) {
            continue;
        }
        const auto nonexc = inc_sed_taxon->getData().nonexcluded_ids;
        assert(nonexc != nullptr);
        auto soln_mrca = find_mrca_without_desids(sampled);
        for (auto t : iter_leaf_n(*soln_mrca)) {
            if (sampled.count(t) == 0 && 0 == nonexc->count(t->getOttId())) {
                unprune_stats.non_monophyletic_inc_sed[inc_sed_taxon] = soln_mrca;
            }
        }
    }

}

// last step in unpruneTaxa
void unprunePreppedInputs(Tree_t & taxonomy,
                          Tree_t & solution,
                          const OttIdSet & incertae_sedis_ids,
                          const map<long, Node_t *> & ott_to_tax,
                          LostTaxonMap & ltm,
                          UnpruneStats & unprune_stats) {
    map<long, Node_t*> ott_to_sol;
    const auto snVec = all_nodes(solution);
    // postorder walk over solution. Every time we find a node assigned to a taxon
    //  we augment the slice of the tree that is rooted at that node (and is the
    //  subtree that is cut at the deepest taxonomic node)
    for (auto nd: snVec){
        if (nd->isTip()) {
            unprune_stats.startNumSolnLeaves += 1;
            if (!nd->hasOttId()) {
                throw OTCError() << "Tip "<< nd->getName() << " in solution lacks an OTT ID";
            }
            if (not ott_to_tax.count(nd->getOttId())) {
                throw OTCError() << "OttId "<< nd->getOttId() << " not in taxonomy!";
            }
        } else {
            if (nd->isOutDegreeOneNode()) {
                unprune_stats.startNumSolnMonotypicInternals += 1;
            } else {
                unprune_stats.startNumSolnForkingInternals += 1;
            }
        }
        auto p = nd->getParent();
        if (nd->hasOttId()){
            ott_to_sol[nd->getOttId()] = nd;
            if (!nd->isTip()) {
                unprune_stats.startNumNamedInternalsInSoln += 1;
                // unpruneTaxaForSubtree will deal with this slice of the tree
                const auto x = unpruneTaxaForSubtree(nd, ott_to_tax, ott_to_sol, ltm);
                unprune_stats.numInternalsAddedToBackbone += std::get<0>(x);
                unprune_stats.numInternalsMergedOnBackbone += std::get<1>(x);
                unprune_stats.numLeavesExpanded += std::get<2>(x);
                if (p) {
                    assert(p == nd->getParent());
                }
            }
            nd->getData().desIds.insert(nd->getOttId());
            if (p) {
                auto & pd = p->getData();
                pd.desIds.insert(nd->getOttId());
            }
        } else {
            if (!nd->isTip()) {
                unprune_stats.numUnnamedInternalsInSoln += 1;
            }
            if (p) {
                const auto & d = nd->getData().desIds;
                p->getData().desIds.insert(d.begin(), d.end());
            }
        }
    }
}

// Note that this function will steal some nodes from `taxonomy`
//  so, on output that will no longer be a valid tree. This should
//  be fine if both taxonomy and solution go out of scope together.
template<typename T>
LostTaxonMap unpruneTaxa(T & taxonomy,
                         T & solution,
                         std::ostream * statsStreamPtr,
                         const OttIdSet & incertae_sedis_ids) {
    LostTaxonMap ltm;
    UnpruneStats unprune_stats;
    unprune_stats.out_degree_many1 = n_internal_out_degree_many(taxonomy);
    // 1. First, index taxonomy by OttId.
    map<long, Node_t *> ott_to_tax;
    OttIdSet monotypicOttIds;
    indexNodesByOttId(taxonomy, ott_to_tax, monotypicOttIds, unprune_stats);
    // filling the nonexcluded IDs (next step) needs the desIds filled for the parts of 
    //    the taxonomy that descend from an incertae sedis taxon.
    // Note that we'll use empty desIds of a taxonomy node as a flag that that taxon is 
    //    not a descendant of an incertae sedis taxon.
    fillDesIdsForIncertaeSedisOnly(ott_to_tax, incertae_sedis_ids, unprune_stats);
    // Now we fill in the unique non-excluded Ids. A non-excluded ID for an internal
    //    node y is the ID of a taxon that is a descendant of some incertae sedis taxon x
    //    where x is child of one of the ancestors of the node y.
    fillNonexcludeIDFields(taxonomy, incertae_sedis_ids);
    // record tips in the sample that are incertae sedis descendants
    findIncertaeSedisInSolution(taxonomy,
                                solution,
                                ott_to_tax,
                                unprune_stats);
    // note which ones can't be in the solution
    findBrokenIncertaeSedisDescendants(taxonomy,
                                       solution,
                                       ott_to_tax,
                                       unprune_stats);
    // Unprune, collecting stats
    unprunePreppedInputs(taxonomy,
                         solution,
                         incertae_sedis_ids,
                         ott_to_tax,
                         ltm,
                         unprune_stats);

    const auto numTaxaMonotypicInternals = monotypicOttIds.size();
    // This is similar to, but different from, the number of non-monotypic nodes reject.
    // That is because the rejected nodes are marked as monotypic if they have no ANCESTRAL children.
    const auto numTaxaRejected = ltm.size();
    size_t numMonotypicTaxaRejected = 0U;
    for (auto monoOtt : monotypicOttIds) {
        if (ltm.count(monoOtt) > 0) {
            numMonotypicTaxaRejected += 1;
        }
    }
    const auto numForkingTaxaRejected = numTaxaRejected - numMonotypicTaxaRejected;
    const auto startNumSolnInternals = unprune_stats.startNumSolnMonotypicInternals + unprune_stats.startNumSolnForkingInternals;
    const auto numTaxaInternals = numTaxaMonotypicInternals + unprune_stats.numTaxaForkingInternals;
    const auto numInternalsAddedByUnpruning = numTaxaInternals - numTaxaRejected - unprune_stats.numInternalsMergedOnBackbone;
    const auto numLeavesAdded = unprune_stats.numTaxaLeaves - unprune_stats.startNumSolnLeaves;
    std::cerr << "Leaves:           solution = " << unprune_stats.startNumSolnLeaves << "   taxonomy = " << unprune_stats.numTaxaLeaves << std::endl;
    std::cerr << "Internal:         solution = " << startNumSolnInternals << "   taxonomy = " << numTaxaInternals << std::endl;
    std::cerr << "Internal splits:  solution = " << unprune_stats.startNumSolnForkingInternals << "   taxonomy = " << unprune_stats.numTaxaForkingInternals << std::endl;
    std::cerr << "Taxa rejected: = " << numTaxaRejected << std::endl;
    std::cerr << "Taxonomy splits: #rejected by phylo inputs = " << numForkingTaxaRejected << std::endl;
    std::cerr << "Solution splits: input # with OTT Ids       = " << unprune_stats.startNumNamedInternalsInSoln << std::endl;
    if (statsStreamPtr) {
        json document;
        json input;
        json output;
        input["num_solution_leaves"] = unprune_stats.startNumSolnLeaves;
        input["num_taxonomy_leaves"] = unprune_stats.numTaxaLeaves;
        input["num_solution_internals"] = startNumSolnInternals;
        input["num_taxonomy_internals"] = numTaxaInternals;
        input["num_solution_splits"] = unprune_stats.startNumSolnForkingInternals;
        input["num_taxonomy_splits"] = unprune_stats.numTaxaForkingInternals;
        input["num_incertae_sedis_taxa"] = incertae_sedis_ids.size();
        output["num_taxa_rejected"] = numTaxaRejected;
        output["num_taxonomy_splits_rejected"] = numForkingTaxaRejected;
        output["num_taxonomy_internals_merged"] = unprune_stats.numInternalsMergedOnBackbone;
        output["num_solution_leaves_expanded"] = unprune_stats.numLeavesExpanded;
        output["num_leaves_added"] = numLeavesAdded;
        document["input"] = input;
        document["output"] = output;
        *statsStreamPtr << document.dump(1) << std::endl;
    }
    return ltm;
}

bool handleLostTaxaJSON(OTCLI&, const string & arg) {
    lostTaxaJSONFilename = arg;
    return true;
}

bool handleStatsJSON(OTCLI&, const string & arg) {
    statsJSONFilename = arg;
    return true;
}

bool handleIncertaeSedis(OTCLI&, const string & arg) {
    incertaeSedisFilename = arg;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-unprune-solution-and-name-unnamed-nodes",
                "Takes a phylogenetic estimate (the first tree file) and a taxonomy (the second tree file).\n"
                "  1. Grafts taxa which were present in the taxonomy, but not included in the solution, onto the.\n"
                "     solution,\n"
                "  2. generates names for the unnamed nodes in the solution\n"
                "  3. writes the augmented solution to standard output as a newick string\n"
                "  4. optionally (if the -j agument is supplied) writes a JSON summary of which higher taxa from\n"
                "     the taxonomy are not included in the solution\n"
                "This program assumes that, if an internal node in the input solution is labelled with a taxonomic\n"
                "identifier, then that taxon is the same definition as the taxon with that identifier in the\n"
                "taxonomic tree. In other words, it does not verify the validity of the internal node names that do exist.\n",
                "-jout/lost-taxa.json solution.tre taxonomy.tre");
    otCLI.addFlag('j',
                  "Produce a JSON file that notes where the descendants of unincluded appear on the tree",
                  handleLostTaxaJSON,
                  true);
    otCLI.addFlag('s',
                  "Produce a JSON file with stats about the inputs and outputs",
                  handleStatsJSON,
                  true);
    otCLI.addFlag('i',
                  "Optional list of IDs of tree in the exemplified taxonomy that are incertae sedis",
                  handleIncertaeSedis,
                  true);
    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};
    if (argc < 2) {
        throw OTCError("No trees provided!");
    }
    if (treeProcessingMain<Tree_t>(otCLI, argc, argv, get, nullptr, 1)) {
        return 1;
    }
    if (trees.size() != 2) {
        std::cerr << "Supplied " << trees.size() << " trees for regrafting, should be 2 trees!";
        return 2;
    }
    std::ofstream lostTaxaJSONStream;
    if (!lostTaxaJSONFilename.empty()) {
        lostTaxaJSONStream.open(lostTaxaJSONFilename);
        if (!lostTaxaJSONStream.good()) {
            std::cerr << "Could not open the path \"" << lostTaxaJSONFilename << "\" to write the JSON file\n";
            return 1;
        }
    }
    std::ofstream statsJSONStream;
    std::ofstream * statsStreamPtr = nullptr;
    if (!statsJSONFilename.empty()) {
        statsJSONStream.open(statsJSONFilename);
        if (!statsJSONStream.good()) {
            std::cerr << "Could not open the path \"" << statsJSONFilename << "\" to write the JSON file\n";
            return 1;
        }
        statsStreamPtr = & statsJSONStream;
    }
    OttIdSet incertae_sedis_ids;
    if (!incertaeSedisFilename.empty()) {
        std::ifstream incert_sed_id_file(incertaeSedisFilename.c_str());
        if (!incert_sed_id_file.good()) {
            std::cerr << "Could not open the incertae sedis ID file: " << incertaeSedisFilename << "\n";
            return 1;
        }
        string line;
        while (std::getline(incert_sed_id_file, line)) {
            char* temp;
            long ott_id = std::strtoul(line.c_str(), &temp, 10);
            if (*temp != '\0' && *temp != '\n') {
                std::cerr << "Expecting just numbers and newlines in incertae sedis file found: " << line << "\n";
                return 1;
            }
            incertae_sedis_ids.insert(ott_id);
        }
    }
    auto & solution = *(trees.at(0));
    auto & taxonomy = *(trees.at(1));
    const auto lostTaxa = unpruneTaxa(taxonomy, solution, statsStreamPtr, incertae_sedis_ids);
    nameUnamedNodes(solution);
    writeTreeAsNewick(std::cout, solution);
    std::cout << std::endl;
    if (!lostTaxaJSONFilename.empty()) {
        writeLostTaxa(lostTaxaJSONStream, lostTaxa);
        lostTaxaJSONStream.close();
    }
    if (statsStreamPtr) {
        statsJSONStream.close();
    }
    return 0;
}
