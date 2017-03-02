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


struct RTNodePartialDesSet;
using Node_t = RootedTreeNode<RTNodePartialDesSet>;

struct RTNodePartialDesSet {
    std::set<long> desIds;
    OttIdSet * nonexcluded_ids = nullptr; // union of IDs that descend from any child of an ancestor marked as incertae sedis
    OttId smallest_child  = 0;
    int tax_level = -1;
    // set for descendants of incertae sedis taxa to indicate which clade hold the taxon
    Node_t * reprInSolnBy = nullptr; 
};

struct RTTreeIncertaeSedisHolder {
    std::list<OttIdSet> ownedIdSets; // used for memory management.
};

using Tree_t = RootedTree<RTNodePartialDesSet, RTTreeIncertaeSedisHolder>;


// Stores pointers for a taxon X which not in the solution.
//   points to the MRCA, and all of the attachment points for 
//   subtrees that are part of this taxon.
// An attachment point is:
//      1. an ancestor of at least one member of the taxon X,
//      2. an ancestor of at least one taxon that is not contained in the taxon X.
//      3. has at least one child (an attached node) that consists entirely
//          of member of taxon X
class LostTaxonDetails {
    using node_vec = std::set<Tree_t::node_type *>;
    using attach2node_vec = std::map<Tree_t::node_type *, node_vec>;
    public:
    LostTaxonDetails(long oid=0, Tree_t::node_type *mrcaNd=nullptr)
        :ottId(oid),
        mrca(mrcaNd) {
    }
    void addAttachmentPoint(Tree_t::node_type *attach, Tree_t::node_type *child) {
        attachNode2AttachedVec[attach].insert(child);
    }
    public: //@TMP should be private.
    long ottId = 0L;
    Tree_t::node_type * mrca = nullptr;
    // for each attachment point, we store which subtrees for this taxon attach there.
    attach2node_vec attachNode2AttachedVec;
    OttIdSet intrudingTaxa;
};
typedef std::map<OttId, LostTaxonDetails> LostTaxonMap;
typedef std::pair<Node_t *, Node_t *> childParPair;

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
    std::set<Node_t *> inc_sed_soln_tips; 
    std::set<Node_t *> inc_sed_internals;
    // maps inc_sed taxa nodes to tips in solution that are descendants of that taxon.
    std::map<Node_t *, std::set<Node_t *> > inc_sed_taxon_to_sampled_tips;
    // maps inc_sed taxon that is not monophyletic to MRCA in solution.
    std::map<Node_t *, Node_t *> non_monophyletic_inc_sed;
    std::map<OttId, OttIdSet> non_monophyletic_inc_sed_to_desIds;
    std::map<OttId, OttIdSet> synonym_map;
    LostTaxonMap lost_taxa;
    OttIdSet monotypic_ott_ids;
    const OttIdSet & incertae_sedis_ids;
    UnpruneStats(const OttIdSet & inc_sed_ids)
        :incertae_sedis_ids(inc_sed_ids) {
    }
};


template<typename T>
void unpruneTaxa(T & taxonomy, T & solution, UnpruneStats & unprune_stats);
void writeLostTaxa(std::ostream & out, const LostTaxonMap & ltm);
template <typename N>
N * find_mrca_without_desids(const std::set<N *> & tip_node_set, std::set<N*> * induced_nodes = nullptr);
template <typename N>
void moveUnsampledChildren(N * taxon, N * solnNode, std::set<N*> *v = nullptr);

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


void writeLostTaxa(ostream & out,
                   const LostTaxonMap & ltm, 
                   const map<OttId, OttIdSet> & synmap) {
    json lostTaxa;
    for (auto ltmEl : ltm) {
        const auto & ottId = ltmEl.first;
        const auto & ltl = ltmEl.second; 
        json lostTaxonLocation;
        lostTaxonLocation["mrca"] = ltl.mrca->get_name();
        const auto attachmentPoints = ltl.attachNode2AttachedVec;
        json apjson;
        for (auto attachPoint: attachmentPoints) {
            const auto & parent = attachPoint.first;
            const auto & childVec = attachPoint.second;
            json childVecJSON = json::array();
            for (auto c : childVec) {
                childVecJSON.push_back(c->get_name());
            }
            apjson[parent->get_name()] = childVecJSON;
        }
        lostTaxonLocation["attachment_points"] = apjson;
        lostTaxa["ott" + std::to_string(ottId)] = lostTaxonLocation;
    }
    json synTaxa;
    for (auto synEl : synmap) {
        const auto & ottId = synEl.first;
        const auto & synset = synEl.second;
        json synj;
        for (auto jsEl : synset) {
            synj.push_back("ott" + std::to_string(jsEl));
        }
        synTaxa["ott" + std::to_string(ottId)] = synj;
    }
    json document;
    document["non_monophyletic_taxa"] = lostTaxa;
    document["taxa_matching_multiple_ott_ids"] = synTaxa;
    out << document.dump(1) << std::endl;
}

// adds ottId to the desIds field of every node from firstNd down to ancAndLast (inclusive)
void addToDesIdsForAnc(long ottId, Node_t *firstNd, Node_t *ancAndLast, const map<long, childParPair> & ott2tax) {
    LOG(DEBUG) << "Adding " << ottId << " to desIds for all taxa from " << firstNd->get_ott_id() << " to " << ancAndLast->get_ott_id();
    assert(firstNd);
    while (true) {
        firstNd->get_data().desIds.insert(ottId);
        //LOG(DEBUG) << "addToDesIdsForAnc hit " << firstNd->get_ott_id();
        //LOG(DEBUG) << "adding " << ottId << " to node " << firstNd->get_name(); dbWriteOttSet(" its desIds", firstNd->get_data().desIds);
        if (firstNd == ancAndLast) {
            return;
        }
        auto fid = firstNd->get_ott_id();
        firstNd = ott2tax.at(fid).second; // parent in unaltered taxonomy.
        // If this triggers, then ancAndLast is not an ancestor of firstNd
        assert(firstNd);
        
    }
}

template <typename N>
vector<N *> higherTaxPreOrderBelowBoundaries(N * root,
                                             const set<N *> & boundaries) {
    set<N *> seen;
    vector<N *> r;
    N * curr = root;
    assert(!curr->get_data().desIds.empty());
    assert(curr->get_data().reprInSolnBy == nullptr);
    stack<N *> toDealWith;
    while (true) {
        LOG(DEBUG) << "higherTaxPreOrderBelowBoundaries visiting " << curr->get_ott_id();
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
                LOG(DEBUG) << "higherTaxPreOrderBelowBoundaries adding " << curr->get_ott_id();
                r.push_back(curr);
                for (auto c : iter_child(*curr)) {
                    auto & cd = c->get_data();
                    if (!cd.desIds.empty() && cd.reprInSolnBy == nullptr) {
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
    nn->setName(taxon->get_name());
    nn->set_ott_id(taxon->get_ott_id());
    taxon->get_data().reprInSolnBy = nn;
    nn->addChild(solnChild);
    nn->get_data().desIds = solnChild->get_data().desIds;
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
void moveUnsampledChildren(N * taxon, N * solnNode, std::set<N*> *v) {
    assert(taxon);
    assert(solnNode);
    const auto children = all_children(taxon);
    for (auto child : children) {
        auto & cd = child->get_data();
        if (cd.desIds.empty() && cd.reprInSolnBy == nullptr) {
            //LOG(DEBUG) << "moveUnsampledChildren moving " << child->get_ott_id();
            child->detachThisNode();
            solnNode->addChild(child);
            if (v) {
                v->insert(child);
            }
        }
    }
}

// Walk tipward from the @mrca until we find subtrees whose
// leaves are purely descendants of @taxon.  Record these
// subtrees as the attachment points in @ltl.
template <typename N>
void registerAttachmentPoints(const N * taxon,
                              N * mrca,
                              LostTaxonDetails & ltl) {
    assert(taxon);
    assert(mrca);
    set<N *> visited;
    stack<N *> toDealWith;
    LOG(DEBUG) << "registerAttachmentPoints for " << taxon->get_ott_id();
    const auto & taxDes = taxon->get_data().desIds;
    assert(taxDes.size() > 1);
    assert(isProperSubset(taxDes, mrca->get_data().desIds));
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
            const auto & solDes = currSolnNd->get_data().desIds;
            assert(solDes != taxDes); // we should not call this if a des of mrca = taxon
            for (auto c : iter_child(*currSolnNd)) {
                const auto & cdesId = c->get_data().desIds;
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


void fixLostTaxonMap(OttId ott_id,
                     const OttIdSet & tip_ids, //@TODO bad name! see note below
                     const map<long, Node_t*> & ott_to_sol,
                     UnpruneStats & unprune_stats) {
    LOG(DEBUG) << "fixLostTaxonMap for " << ott_id;
    dbWriteOttSet("  tip_ids", tip_ids);
    std::set<Node_t *> tips;
    for (auto i : tip_ids) {
        auto x = ott_to_sol.find(i);
        // "tip_ids" currently is all desIds of the incertae sedis taxon  - includes higher taxa which may not be in soln.
        if (x != ott_to_sol.end()) {
            tips.insert(x->second);
        }
    }
    std::set<Node_t *> visited;
    Node_t * mrca = find_mrca_without_desids(tips, &visited);
    std::set<Node_t *> no_interloper;
    for (auto n: visited) {
        if (tips.count(n) == 0) {
            bool has_interloper = false;
            for (auto c : iter_child(*n)) {
                if (visited.count(c) == 0) {
                    has_interloper = true;
                    break;
                }
            }
            if (!has_interloper) {
                no_interloper.insert(n);
            }
        }
    }
    set<Node_t *> at_set;
    for (auto n : tips) {
        if (n == mrca) {
            at_set.insert(mrca);
            continue;
        }
        auto p = n->getParent();
        while (p != mrca && no_interloper.count(p) > 0) {
            p = p->getParent();
        }
        at_set.insert(p);
    }
    LostTaxonDetails ltl(ott_id, mrca);
    unprune_stats.lost_taxa[ott_id] = ltl;
    auto & attachment = unprune_stats.lost_taxa[ott_id].attachNode2AttachedVec;
    attachment.clear();
    for (auto n: at_set) {
        auto & att_vec = attachment[n];
        for (auto c : iter_child(*n)) {
            if (tips.count(c) >0 || no_interloper.count(c) > 0) {
                att_vec.insert(c);
            }
        }
    }
}

template <typename N>
void addChildrenOfNonMonophyleticTaxon(N * taxon,
                                       N * mrca,
                                       LostTaxonMap & ltm,
                                       bool registerInLTM) {
    assert(taxon);
    assert(mrca);
    assert(taxon->has_ott_id());
    const auto ottId = taxon->get_ott_id();
    LostTaxonDetails * ltl = nullptr;
    if (registerInLTM) {
        if (ltm.count(ottId) > 0) {
            LOG(DEBUG) << " LTM for " << ottId << " again. this can happen if it is incertae sedis.";
        }
        ltm[ottId] = LostTaxonDetails(ottId, mrca);
        ltl = &(ltm[ottId]);
        registerAttachmentPoints(taxon, mrca, *ltl);
    }
    std::set<N *> v;
    moveUnsampledChildren(taxon, mrca, &v);
    if (registerInLTM && !v.empty()) {
        auto & s = ltl->attachNode2AttachedVec[mrca];
        s.insert(v.begin(), v.end());
    }
}


template <typename N>
N * special_bisect_with_new_child(N * taxNode, N * currSolnNd, bool moveUnsamp) {
    N * nt = bisect_branch_with_new_child(currSolnNd);
    nt->setName(taxNode->get_name());
    nt->set_ott_id(taxNode->get_ott_id());
    nt->get_data().desIds = currSolnNd->get_data().desIds;
    if (moveUnsamp) {
        moveUnsampledChildren(taxNode, nt);
    }
    return nt;
}

// Looks through the children of `currSolnNd` If one has all the desIds that are
//    flagged by taxon_node, it returns that node. Otherwise nullptr.
// Used for starting from an ancestor of the taxa in taxon_node and moving one
//    step closer to the MRCA of those taxa
Node_t * find_single_child_with_all_marked_taxa(Node_t * currSolnNd, const Node_t * taxon_node) {
    dbWriteOttSet(" find_single_child_with_all_marked_taxa taxon desIds", taxon_node->get_data().desIds);
    dbWriteOttSet(" find_single_child_with_all_marked_taxa currSolnNd desIds", currSolnNd->get_data().desIds);
    const auto & taxon_data = taxon_node->get_data();
    const auto & taxon_des = taxon_data.desIds;
    Node_t * nextSolnNd = nullptr;
    for (auto c : iter_child(*currSolnNd)) {
        const auto & cdesId = c->get_data().desIds;
        if (!areDisjoint(cdesId, taxon_des)) {
            //LOG(DEBUG) << "not disjoint nextSolnNd=" << (long) nextSolnNd;
            if (nextSolnNd == nullptr) {
                nextSolnNd = c;
            } else {
                return nullptr;
            }
        }
    }
    return nextSolnNd;
}

// returns the number of nodes on the backbone that are merged with this higher taxon
//      because the taxon was compatible (the same as) the node - will be either 0 or 1.
template <typename N>
size_t incorporateHigherTaxonNode(N* higherTaxonNd,
                                  N* rootSolnNd,
                                  set<N *> & nodesAddedForTaxa,
                                  UnpruneStats & unprune_stats,
                                  bool recordInLTM,
                                  const map<long, childParPair> & ott2tax) {
    LostTaxonMap & ltm = unprune_stats.lost_taxa;
    assert(higherTaxonNd);
    assert(rootSolnNd);
    LOG(DEBUG) << "higherTaxonNd->name = " << higherTaxonNd->get_name();
    const auto & taxData = higherTaxonNd->get_data();
    const auto & taxDes = taxData.desIds;
    dbWriteOttSet(" higherTaxonNd desIds", taxDes);
    assert(taxDes.size() > 0);
    assert(taxData.nonexcluded_ids != nullptr);
    const auto & nonexcluded_ids = *taxData.nonexcluded_ids;
    N * currSolnNd = rootSolnNd;
    N * tipmostMRCA = nullptr;
    // the root of the solution, must have all of the constituent taxa
    assert(isSubset(taxDes, rootSolnNd->get_data().desIds));
    // here we walk tipward to find the MRCA of the taxa in taxDes
    while (true) {
        N * nextSolnNd = find_single_child_with_all_marked_taxa(currSolnNd, higherTaxonNd);
        if (nextSolnNd == nullptr) {
            break;
        }
        currSolnNd = nextSolnNd;
    }
    if (unprune_stats.non_monophyletic_inc_sed.count(higherTaxonNd) > 0) {
        addChildrenOfNonMonophyleticTaxon(higherTaxonNd, currSolnNd, ltm, recordInLTM);
    }
    assert(currSolnNd != nullptr);
    const auto & solDes = currSolnNd->get_data().desIds;
    assert(isSubset(taxDes, solDes));
    auto extra = set_difference_as_set(solDes, taxDes);
    bool hasExcludedExtras = false;
    for (auto e: extra) {
        if (nonexcluded_ids.count(e) == 0) {
            hasExcludedExtras = true;
            break;
        }
    }
    if (hasExcludedExtras) {
        // the MRCA of this taxon in the solution contains nodes that are in the taxon's exclude
        //    set. So this taxon in broken in this solution.
        LOG(DEBUG) << "taxon " << higherTaxonNd->get_ott_id() << " is not monophyletic. Adding other children...";
        addChildrenOfNonMonophyleticTaxon(higherTaxonNd, currSolnNd, ltm, recordInLTM);
        return 0U;
    }
    // This solution node has the same descendant set of this taxon.
    // It could either be the only child of the taxon (if the taxon) is 
    //    monotypic, or it could be labelled with this OTT ID.
    N * nt = nullptr;
    // the easiest case is if the solution node has the same taxon ID
    //    just move the unsampled children.
    if (currSolnNd->has_ott_id() && currSolnNd->get_ott_id() == higherTaxonNd->get_ott_id()) {
        LOG(DEBUG) << "same ID call moveUnsampledChildren";
        moveUnsampledChildren(higherTaxonNd, currSolnNd);
        return 1U;
    }
    // The monotypic case is also easy... add a monotypic node to the solution.
    if (higherTaxonNd->isOutDegreeOneNode()) {
        // careful to keep the root of the subproblem the same node structure, which
        //    may require swapping data to make the root the new deeper monotypic taxon....
        if (currSolnNd == rootSolnNd) {
            LOG(DEBUG) << "adding Monotypic child";
            nt = special_bisect_with_new_child(higherTaxonNd, currSolnNd, false);
        } else {
            LOG(DEBUG) << "introduceMonotypicParent";
            nt = introduceMonotypicParent(higherTaxonNd, currSolnNd);
        }
    } else {
        // In the world of incertae sedis, we need to worry about the case of the node
        //    already having a label, but not being a descendant of higher taxon.
        // We only want to add a new node to the solution if the solution node that
        //    we have found is associated with a descendant taxon of
        bool add_new_node = false;
        if (currSolnNd->has_ott_id()) {
            // we have to check the taxonomy tree to find out...
            auto potential_des = ott2tax.at(currSolnNd->get_ott_id()).first;
            add_new_node = isAncDecPair(higherTaxonNd, potential_des);
            if ((!add_new_node) && rootSolnNd == currSolnNd) {
                add_new_node = isAncDecPair(potential_des, higherTaxonNd);
            }
        }
        if (add_new_node) {
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
            if (currSolnNd->has_ott_id()) {
                unprune_stats.synonym_map[currSolnNd->get_ott_id()].insert(higherTaxonNd->get_ott_id());
            } else {
                currSolnNd->setName(higherTaxonNd->get_name());
                currSolnNd->set_ott_id(higherTaxonNd->get_ott_id());
            }
        }
    }
    if (nt != nullptr) {
        nodesAddedForTaxa.insert(nt);
        return 0U;
    }
    return 1U;
}

typedef pair<Node_t *, Node_t *> TaxSolnNdPair;

void registerCoveringOfIncSed(OttId effectiveTipOttId, 
                              Node_t * effTipTaxonNd,
                              const map<long, Node_t *> & ott2soln, 
                              const UnpruneStats & unprune_stats,
                              std::map<Node_t *, std::set<TaxSolnNdPair> > & curr_slice_inc_sed_map,
                              const map<long, childParPair> & ott2tax) {
    //LOG(DEBUG) << "registerCoveringOfIncSed for " << effectiveTipOttId;
    const auto & inc_sed_internals = unprune_stats.inc_sed_internals;
    auto anc = ott2tax.at(effectiveTipOttId).second;
    effTipTaxonNd->get_data().desIds.insert(effectiveTipOttId);
    auto effTipSolnNd = ott2soln.at(effectiveTipOttId);
    if (inc_sed_internals.find(effTipTaxonNd) != inc_sed_internals.end()) {
        curr_slice_inc_sed_map[effTipTaxonNd].insert(TaxSolnNdPair(effTipTaxonNd, effTipSolnNd));
    }
    while (anc &&  inc_sed_internals.find(anc) != inc_sed_internals.end()) {
        //LOG(DEBUG) << "registerCoveringOfIncSed added " << effectiveTipOttId << " to " << anc->get_ott_id();
        anc->get_data().desIds.insert(effectiveTipOttId);
        curr_slice_inc_sed_map[anc].insert(TaxSolnNdPair(effTipTaxonNd, effTipSolnNd));
        anc = ott2tax.at(anc->get_ott_id()).second;
    }
}

// Moves unsampled taxa from the taxanomy (by finding them in `ott2tax`) to a slice of the
//  solution tree that is rooted at `rootSolnNd`.
// Uses the `rootSolnNd->get_data().desIds` to denote the set of leaves of this slice of the
//  tree. On exit this will be the ID of this clade and the IDs of any uncoalesced incertae sedis 
//    taxa that are in this slice.
// If the taxon conflicts with the solution, then the unsampled taxa are attached at the
//  MRCA of the taxon in the solution, and an entry is added mapping the taxon's OTT Id
//  to a LostTaxonDetails object that explains where the elements of the taxon are located.
// Returns the set of IDs that should  
template <typename N>
void unpruneTaxaForSubtree(N *rootSolnNd,
                           const map<long, childParPair> & ott2tax, 
                           map<long, N*> & ott2soln, 
                           UnpruneStats & unprune_stats) {
    assert(rootSolnNd);
    assert(rootSolnNd->has_ott_id());
    const auto ottId = rootSolnNd->get_ott_id();
    LOG(DEBUG) << "unpruneTaxaForSubtree for " << rootSolnNd->get_name();
    assert(ott2soln.at(ottId) == rootSolnNd);
    N * rootTaxonNd = ott2tax.at(ottId).first;
    assert(rootTaxonNd != nullptr);
    assert(!rootTaxonNd->isTip());
    const OttIdSet * rootNonExclude = rootTaxonNd->get_data().nonexcluded_ids;
    assert(rootNonExclude != nullptr);
    // this will have the IDs for all of the include taxa for this slice of the tree
    // note that some may be "higher" taxa.
    auto & solnDesIds = rootSolnNd->get_data().desIds;
    assert(!solnDesIds.empty());
    // desIds fields of the solution tree are filled in by the caller before this function.
    //  here we fill in those fields for the taxonomy nodes.
    //Here we add the sampled IDs to desId fields ofthe relevant taxa. This is potentially
    //  confusing because get_data().desIds will hold only the IDs of taxa that are present in the soln
    //  tree (even when we are talking about the taxonomy node's desIds)
    // While we are walking through the "leaf" Ids for this tree slice, we'll also
    //  collect the leaf sets
    const auto & inc_sed_map = unprune_stats.inc_sed_taxon_to_sampled_tips;
    const auto & inc_sed_internals = unprune_stats.inc_sed_internals;
    set<N *> solnLeaves;
    set<N *> taxaLeaves;
    std::map<Node_t *, std::set<TaxSolnNdPair> > curr_slice_inc_sed_map;
    set<Node_t *> inc_sed_that_are_proper_children;
    dbWriteOttSet(" effective tips: ", solnDesIds);
    dbWriteOttSet(" rootNonExclude: ", *rootNonExclude);
    list<Node_t *> effTipTaxa;
    typedef pair<OttId, OttId> carriedIDRealId;
    list<carriedIDRealId > queuedForFlagging;
    for (auto effectiveTipOttId : solnDesIds) {
        auto effTipTaxonNd = ott2tax.at(effectiveTipOttId).first;
        effTipTaxa.push_back(effTipTaxonNd);
        auto is_map_it = inc_sed_map.find(effTipTaxonNd);
        // see if this "tip" is an incertae sedis taxon.
        bool is_inc_sed = is_map_it != inc_sed_map.end();
        bool is_tip_inc_sed = false;
        if (!is_inc_sed && effTipTaxonNd->isTip()) {
            auto corres_soln_tip = ott2soln.at(effectiveTipOttId);
            if (unprune_stats.inc_sed_soln_tips.count(corres_soln_tip) > 0) {
                is_inc_sed = true;
                is_tip_inc_sed = true;
            }
        }
        bool is_proper_des = false;
        if (is_inc_sed) {
            if (rootNonExclude->count(effectiveTipOttId) > 0) {
                LOG(DEBUG) << "checking effective tip " << effectiveTipOttId << " lower  + inc_sed";
                // this "tip" taxon is nonexcluded incertae sedis taxon from deeper in tree...
                if (!is_tip_inc_sed && is_map_it->second.size() > 0) {
                    throw OTCError() << "Incertae sedis taxon " << effectiveTipOttId << " is a tip for slice " << ottId << ", but is still in the inc_sed_taxon_to_sampled_tips map!";
                }
                registerCoveringOfIncSed(effectiveTipOttId, 
                                         effTipTaxonNd,
                                         ott2soln, 
                                         unprune_stats,
                                         curr_slice_inc_sed_map,
                                         ott2tax);
                if (effTipTaxonNd->isTip()) {
                    taxaLeaves.insert(effTipTaxonNd);
                }
            } else {
                LOG(DEBUG) << "checking effective tip " << effectiveTipOttId << " proper +  inc_sed";
                is_proper_des = true;
            }
        } else {
            LOG(DEBUG) << "checking effective tip " << effectiveTipOttId << " not inc_sed";
            is_proper_des = true;
        }
        if (is_proper_des) {
            if (is_inc_sed && !is_tip_inc_sed) {
                inc_sed_that_are_proper_children.insert(effTipTaxonNd);
            }
            if (effTipTaxonNd->get_data().reprInSolnBy != nullptr) {
                assert(effTipTaxonNd->get_data().reprInSolnBy->has_ott_id());
                queuedForFlagging.push_back(carriedIDRealId(effectiveTipOttId, effTipTaxonNd->get_data().reprInSolnBy->get_ott_id()));
                LOG(DEBUG) << "checking effective tip " << effectiveTipOttId << " threaded through " << effTipTaxonNd->get_data().reprInSolnBy->get_ott_id();
            } else {
                queuedForFlagging.push_back(carriedIDRealId(effectiveTipOttId, effectiveTipOttId));
            }
        }    
    }
    bool root_taxon_is_inc_sed = inc_sed_map.find(rootTaxonNd) != inc_sed_map.end();
    if (root_taxon_is_inc_sed) {
        inc_sed_that_are_proper_children.insert(rootTaxonNd);
    }
    // now we flag all of the desIds for the proper
    for (auto cr_pair : queuedForFlagging) {
        auto realTipOttId = cr_pair.second;
        auto effTipTaxonNd = ott2tax.at(realTipOttId).first;
        effTipTaxonNd->get_data().desIds.clear();
    }
    for (auto cr_pair : queuedForFlagging) {
        auto carriedTipId = cr_pair.first;
        auto threadedThroughOttId = cr_pair.second;
        Node_t * effTipTaxonNd = ott2tax.at(carriedTipId).first;
        addToDesIdsForAnc(carriedTipId, effTipTaxonNd, rootTaxonNd, ott2tax);
        solnLeaves.insert(ott2soln.at(threadedThroughOttId));
        taxaLeaves.insert(effTipTaxonNd);
    }   
    // need to process inc. sed. splits in reverse order by level, to mimic the postorder sweep
    map<int, list<Node_t *> > inc_sed_to_deal_with_by_level;
    set<Node_t *> inc_sed_mapping_deeper;

    // DEBUGGING
    OttIdSet db_inc_sed_ott_ids;
    for (auto scism_it: curr_slice_inc_sed_map) {
        auto & inc_sed_taxon = scism_it.first;
        db_inc_sed_ott_ids.insert(inc_sed_taxon->get_ott_id());
    }
    dbWriteOttSet("Incertae sedis taxa for the current slice...", db_inc_sed_ott_ids);
    // End DB out

    // Now we can use curr_slice_inc_sed_map to figure out which inc. sedis. taxa will be dealt with in this slice;
    for (auto scism_it: curr_slice_inc_sed_map) {
        auto & inc_sed_taxon = scism_it.first;
        auto & included_leaf_pair_set = scism_it.second;
        auto gsit = inc_sed_map.find(inc_sed_taxon);
        if (gsit == inc_sed_map.end()) {
            throw OTCError() << "Incertae sedis taxon " << inc_sed_taxon->get_ott_id() << " included in slice " << ottId << ", but was not found in the inc_sed_taxon_to_sampled_tips map!";
        }
        auto & full_soln_set = gsit->second;
        if (full_soln_set.size() == 0) {
            continue; // we don't need to worry about incertae sedis taxa with no descendants stored. these have been dealt with
        }
        LOG(DEBUG) << "inc. sed. " << inc_sed_taxon->get_ott_id() << " in this slice with " << included_leaf_pair_set.size() << " of " << gsit->second.size() << " sampled tips";
        if (full_soln_set.size() == included_leaf_pair_set.size()) {
            // all of the sampled member of this incertae sedis taxon are sampled in this slice...
            int tax_level = inc_sed_taxon->get_data().tax_level;
            inc_sed_to_deal_with_by_level[tax_level].push_back(inc_sed_taxon);
        } else {
            // on exit we need to deal up date outgoing edges...
            inc_sed_mapping_deeper.insert(inc_sed_taxon);
        }
    }

    size_t numExpanded = 0;
    // If a tip of the solution is a higher taxon, then we should
    //  graft on the other tips here...
    for (auto l : solnLeaves) {
        if (l->isTip()) {
            assert(l->has_ott_id());
            auto leafOttId = l->get_ott_id();
            auto taxonForLeaf = ott2tax.at(leafOttId).first;
            if (!taxonForLeaf->isTip()) {
                numExpanded += 1;
                auto cvec = all_children(taxonForLeaf);
                for (auto c : cvec) {
                    if (c->get_data().reprInSolnBy == nullptr) {
                        taxonForLeaf->removeChild(c);
                        c->get_data().reprInSolnBy = l;
                        l->addChild(c);
                    }
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
    auto postOrderInTaxNd = higherTaxPreOrderBelowBoundaries(rootTaxonNd, taxaLeaves);
    std::reverse(postOrderInTaxNd.begin(), postOrderInTaxNd.end());
    // Here we will add the relevant higher (non-leaf) taxa from the taxonomy to the 
    //  solution. It is crucial that we do this in postorder because if we have a series
    //  of nodes to introduce along a branch, the incorporateHigherTaxonNode function
    //  will add them below the attachment node (so postorder will assure that they 
    //  are correctly added in the tip->root orientation).
    set<N *> nodesAddedForTaxa;
    for (auto higherTaxonNd : postOrderInTaxNd) {
        LOG(DEBUG) << "dealing with " << higherTaxonNd->get_ott_id();
        if (higherTaxonNd == rootTaxonNd) {
            moveUnsampledChildren(higherTaxonNd, rootSolnNd);
            numMerged += 1 ;
        } else {
            bool recordInLTM = inc_sed_that_are_proper_children.find(higherTaxonNd) == inc_sed_that_are_proper_children.end();
            numMerged += incorporateHigherTaxonNode(higherTaxonNd, rootSolnNd, nodesAddedForTaxa, unprune_stats, recordInLTM, ott2tax);
        }
        LOG(DEBUG) << "reprInSolnBy marked for " << higherTaxonNd->get_ott_id();
        higherTaxonNd->get_data().reprInSolnBy = rootSolnNd;
    }
    for (auto is_it = inc_sed_to_deal_with_by_level.rbegin(); is_it != inc_sed_to_deal_with_by_level.rend(); ++is_it) {
        for (auto is_taxon_ptr : is_it->second) {
            LOG(DEBUG) << "dealing with Incertae sedis " << is_taxon_ptr->get_ott_id() << " from level " << is_it->first;
            if (is_taxon_ptr == rootTaxonNd) {
                moveUnsampledChildren(is_taxon_ptr, rootSolnNd);
                numMerged += 1 ;
            } else {
                numMerged += incorporateHigherTaxonNode(is_taxon_ptr, rootSolnNd, nodesAddedForTaxa, unprune_stats, false, ott2tax);
            }
            LOG(DEBUG) << "reprInSolnBy marked for " << is_taxon_ptr->get_ott_id();
            is_taxon_ptr->get_data().reprInSolnBy = rootSolnNd;
        }
    }
    // for the sake of a low memory footprint, here we 
    //  clear out the desIds fields for this slice of the tree.
    for (auto n : taxaLeaves) {
        n->get_data().reprInSolnBy = rootSolnNd;
        n->get_data().desIds.clear();
    }
    for (auto n : postOrderInTaxNd) {
        n->get_data().desIds.clear();
    }
    set<N *> visited;
    for (auto n : solnLeaves) {
        auto c = n;
        while (visited.count(c) == 0 && c != rootSolnNd) {
            c->get_data().desIds.clear();
            if (c != n) {
                visited.insert(c);
            }
            c = c->getParent();
        }
    }
    for (auto scism_it: curr_slice_inc_sed_map) {
        auto & inc_sed_taxon = scism_it.first;
        inc_sed_taxon->get_data().desIds.clear();
    }
    rootSolnNd->get_data().desIds.clear();
    rootSolnNd->get_data().desIds.insert(ottId);
    // here we need to update the altered inc_sed taxon maps...
    // All the ones finished in this slice can be deleted from maps...
    for (auto is_taxon_ptr : inc_sed_that_are_proper_children) {
        LOG(DEBUG) << "removing des " << is_taxon_ptr->get_ott_id() << " from inc_sed_taxon_to_sampled_tips";
        unprune_stats.inc_sed_taxon_to_sampled_tips.at(is_taxon_ptr).clear();
        if (is_taxon_ptr->getParent()) {
            is_taxon_ptr->detachThisNode();
        }
    }
    for (auto is_it = inc_sed_to_deal_with_by_level.rbegin(); is_it != inc_sed_to_deal_with_by_level.rend(); ++is_it) {
        for (auto is_taxon_ptr : is_it->second) {
            LOG(DEBUG) << "removing " << is_taxon_ptr->get_ott_id() << " from inc_sed_taxon_to_sampled_tips";
            unprune_stats.inc_sed_taxon_to_sampled_tips.at(is_taxon_ptr).clear();
            if (is_taxon_ptr->getParent()) {
                is_taxon_ptr->detachThisNode();
            }
        }
    }
    // all of the ones with MRCA deeper need to have the rootSolnNd registered as the "sampled tip" (even though it really isn't a tip)
    for (auto is_taxon_ptr : inc_sed_mapping_deeper) {
        auto & samp_tip_set = unprune_stats.inc_sed_taxon_to_sampled_tips.at(is_taxon_ptr);
        LOG(DEBUG) << "altering " << is_taxon_ptr->get_ott_id() << " from inc_sed_taxon_to_sampled_tips from " << samp_tip_set.size();
        auto & stsp = curr_slice_inc_sed_map[is_taxon_ptr];
        for (TaxSolnNdPair tsp : stsp) {
            // for the deeper parts of the tree the tips in this slice need to be annotated
            //    to reflect the fact that they have already been included in a taxonomic node
            Node_t * spike_nd = tsp.first;
            rootSolnNd->get_data().desIds.insert(spike_nd->get_ott_id());
            //flagIncSedAs
            while (true) {
                if (inc_sed_mapping_deeper.count(spike_nd) > 0) {
                    break;
                }
                LOG(DEBUG) << "reprInSolnBy marked for " << spike_nd->get_ott_id();
                spike_nd->get_data().reprInSolnBy = rootTaxonNd;
                spike_nd = spike_nd->getParent();
                assert(spike_nd != nullptr);
            }
            samp_tip_set.erase(tsp.second);
        }
        samp_tip_set.insert(rootSolnNd);
        LOG(DEBUG) << "    to " << samp_tip_set.size();
    }
    unprune_stats.numInternalsAddedToBackbone += nodesAddedForTaxa.size();
    unprune_stats.numInternalsMergedOnBackbone += numMerged;
    unprune_stats.numLeavesExpanded += numExpanded;
}

void fillNonexcludeIDFields(Tree_t & taxonomy, const OttIdSet & incertae_sedis_ids) {
    // Assuming that incertae sedis taxa are fairly rare, there will be lots of nodes that 
    //    share the same set of non-excluded taxa. Thus we have an owning store and lots
    //    of pointers to them.
    auto & owningList = taxonomy.get_data().ownedIdSets;
    auto taxonomy_root_ptr = taxonomy.getRoot();
    owningList.push_back(OttIdSet());
    taxonomy_root_ptr->get_data().nonexcluded_ids = &(*owningList.rbegin());
    taxonomy_root_ptr->get_data().tax_level = 0;
    for (auto nd: iter_pre(taxonomy)) {
        if (nd->isTip()) {
            continue;
        }
        auto & ndd = nd->get_data();
        auto * curr_nonexc = ndd.nonexcluded_ids;
        assert(curr_nonexc != nullptr);
        set<Node_t *> inc_sed_children;
        set<Node_t *> normal_children;
        int child_level = 1 + ndd.tax_level;
        for (auto c : iter_child(*nd)) {
            c->get_data().tax_level = child_level;
            if (incertae_sedis_ids.count(c->get_ott_id())) {
                inc_sed_children.insert(c);
            } else {
                normal_children.insert(c);
            }
        }
        if (inc_sed_children.empty()) {
            // all children can share same nonexcluded set as curr nd;
            for (auto c: normal_children) {
                c->get_data().nonexcluded_ids = curr_nonexc;
            }
        } else {
            if (inc_sed_children.size() == 1) {
                // if there is only one, the incertae sedis child
                //     actually just has the same non-excluded as parent
                auto inc_sed_child = *inc_sed_children.begin();
                auto & inc_sed_child_data = inc_sed_child->get_data();
                inc_sed_child_data.nonexcluded_ids = curr_nonexc;
            } else {
                // each incertae sedis child will have its own non-excluded set
                for (auto isc: inc_sed_children) {
                    owningList.push_back(*curr_nonexc);
                    isc->get_data().nonexcluded_ids = &(*owningList.rbegin());
                }
                for (auto isc: inc_sed_children) {
                    const auto & isc_di = isc->get_data().desIds;
                    for (auto oisc: inc_sed_children) {
                        if (oisc == isc) {
                            continue;
                        }
                        oisc->get_data().nonexcluded_ids->insert(isc_di.begin(), isc_di.end());
                    }
                }
            }
            // all "normal" children share non-excluded, but it is larger than current node's
            if (!normal_children.empty()) {
                owningList.push_back(*curr_nonexc);
                auto normal_child_non_exc = &(*owningList.rbegin());
                for (auto isc: inc_sed_children) {
                    auto isc_di = isc->get_data().desIds;
                    normal_child_non_exc->insert(isc_di.begin(), isc_di.end());
                }
                for (auto c: normal_children) {
                    c->get_data().nonexcluded_ids = normal_child_non_exc;
                }
            }
        }
    }
}

void fillDesIdsForIncertaeSedisOnly(const map<long, childParPair> & ott_to_tax,
                                    const OttIdSet & incertae_sedis_ids, 
                                    UnpruneStats & unprune_stats) {
    // As a part of implementing the incertae sedis logic, we also fill the desIds of the taxonomy.
    // This info is used to create the "non-excluded" sets, but it is expensive to do on the whole 
    //    taxonomy, so we just do it for subtrees rooted at an incertae sedis taxon...
    for (auto ist_id: incertae_sedis_ids) {
        //LOG(DEBUG) << ""
        if (ott_to_tax.count(ist_id) == 0) {
            throw OTCError() << "Incertae sedis ID (" << ist_id << ") not in taxonomy.";
        }
        LOG(DEBUG) << "Filling desIds for incertae sedis " << ist_id;
        auto taxon = ott_to_tax.at(ist_id).first;
        for (auto nd: iter_post_n(*taxon)) {
            auto & di_ref = nd->get_data().desIds;
            // empty check is just an optimization to avoid recalculating the results for a node.
            if (di_ref.empty()) {
                di_ref.insert(nd->get_ott_id());
                unprune_stats.inc_sed_internals.insert(nd);
                for (auto c : iter_child(*nd)) {
                    const auto & cdi_ref = c->get_data().desIds;
                    di_ref.insert(cdi_ref.begin(), cdi_ref.end());
                }
            }
        }
        dbWriteOttSet("the desIds set is ", taxon->get_data().desIds);
    }
}


void indexNodesByOttId(Tree_t & taxonomy,
                       map<long, childParPair> & ott_to_tax,
                       UnpruneStats & unprune_stats) {
    OttIdSet & moi = unprune_stats.monotypic_ott_ids;
    for (auto nd: iter_post(taxonomy)){
        if (nd->has_ott_id()) {
            const OttId taxon_ott_id = nd->get_ott_id();
            ott_to_tax[taxon_ott_id] = childParPair(nd, nd->getParent());
            if (nd->isOutDegreeOneNode()) {
                moi.insert(nd->get_ott_id());
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
                                 const map<long, childParPair> & ott_to_tax,
                                 UnpruneStats & unprune_stats) {
    auto & to_tips_map = unprune_stats.inc_sed_taxon_to_sampled_tips;
    for (auto snd: iter_leaf(solution)) {
        if (!snd->has_ott_id()) {
            throw OTCError() << "Tip "<< snd->get_name() << " in solution lacks an OTT ID";
        }
        auto ott_id = snd->get_ott_id();
        auto taxon_node = ott_to_tax.at(ott_id).first;
        if (ott_id == 5144555) {
            dbWriteOttSet("target desIds", taxon_node->get_data().desIds);
        }
        if (taxon_node->get_data().desIds.empty()) {
            continue;    
        }
        unprune_stats.inc_sed_soln_tips.insert(snd);
        auto anc = taxon_node->getParent();
        to_tips_map[taxon_node].insert(snd);
        while (anc && !(anc->get_data().desIds.empty())) {
            to_tips_map[anc].insert(snd);
            anc = anc->getParent();
        }
    }
    // db
    for (auto x : to_tips_map) {
        LOG(DEBUG) << "findIncertaeSedisInSolution recorded entry in inc_sed_taxon_to_sampled_tips for " << x.first->get_ott_id();    
    }
}

template <typename N>
N * find_mrca_without_desids(const set<N *> & tip_node_set, set<N*> * induced_nodes) {
    if (tip_node_set.size() < 2)  {
        if (tip_node_set.empty()) {
            return nullptr;
        }
        N * r = *tip_node_set.begin();
        if (induced_nodes) {
            induced_nodes->insert(r);
        }
        return r;
    }
    auto nit = tip_node_set.begin();
    N * curr = *nit++;
    deque<N *> root_to_mrca;
    set<N *> local_induced_nodes;
    set<N *> * seen_nodes = (induced_nodes == nullptr ? &local_induced_nodes : induced_nodes);
    set<N *> dequed_nodes;
    while (curr != nullptr) {
        // if (curr->has_ott_id())
        //     LOG(DEBUG) << "mrca pushing " << curr->get_ott_id();
        // else
        //     LOG(DEBUG) << "mrca pushing address " << curr;
        root_to_mrca.push_front(curr);
        seen_nodes->insert(curr);
        dequed_nodes.insert(curr);
        curr = curr->getParent();
    }
    for (; nit != tip_node_set.end(); ++nit) {
        curr = *nit;
        // if (curr->has_ott_id())
        //     LOG(DEBUG) << "mrca checking " << curr->get_ott_id();
        // else
        //     LOG(DEBUG) << "mrca checking address " << curr;
        while (curr != nullptr) {
            if (seen_nodes->count(curr) != 0) {
                if (dequed_nodes.count(curr) != 0) {
                    while (root_to_mrca.back() != curr) {
                        auto td = root_to_mrca.back();
                        // if (td->has_ott_id())
                        //     LOG(DEBUG) << "mrca popping " << td->get_ott_id();
                        // else
                        //     LOG(DEBUG) << "mrca popping address " << td;
                        
                        root_to_mrca.pop_back();
                        dequed_nodes.erase(td);
                        if (root_to_mrca.size() == 1) {
                            return root_to_mrca.back();
                        }
                    }
                }
                break;
            } else {
                seen_nodes->insert(curr);
            }
            curr = curr->getParent();
        }
    }
    N * r = root_to_mrca.back();
    if (induced_nodes != nullptr) {
        N * p = r->getParent();
        while (p != nullptr) {
            if (induced_nodes->count(p) > 0) {
                induced_nodes->erase(p);
            } else {
                break;
            }
            p = p->getParent();
        }
    }
    return r;
}

// Taxa that are descendants of an incertae sedis taxon, will demand some special handling..
// so we find out which ones are not found in the solution before calling unprunePreppedInputs
// results are stored in unprune_stats. This relies on having been filled...
void findBrokenIncertaeSedisDescendants(Tree_t & taxonomy,
                                        Tree_t & solution,
                                        const map<long, childParPair> & ott_to_tax,
                                        UnpruneStats & unprune_stats) {
    const auto & to_tips_map = unprune_stats.inc_sed_taxon_to_sampled_tips;
    LOG(DEBUG) << " # incertae_sedis is " << unprune_stats.inc_sed_taxon_to_sampled_tips.size();
    for (auto ttm_it: to_tips_map) {
        auto inc_sed_taxon = ttm_it.first;
        LOG(DEBUG) << " incertae_sedis taxon: " << inc_sed_taxon->get_ott_id();
        const auto & sampled = ttm_it.second;
        if (sampled.size() < 2) {
            continue;
        }
        const auto nonexc = inc_sed_taxon->get_data().nonexcluded_ids;
        const auto taxon_des = inc_sed_taxon->get_data().desIds;
        assert(nonexc != nullptr);
        auto soln_mrca = find_mrca_without_desids(sampled);
        for (auto t : iter_leaf_n(*soln_mrca)) {
            LOG(DEBUG) << " visiting leaf " << t->get_ott_id() << " to see if it breaks " << inc_sed_taxon->get_ott_id();
            if (sampled.count(t) == 0 && 0 == nonexc->count(t->get_ott_id())) {
                LOG(DEBUG) << "Inc. sed. " << inc_sed_taxon->get_ott_id() << " is not monophyletic.";
                unprune_stats.non_monophyletic_inc_sed[inc_sed_taxon] = soln_mrca;
                unprune_stats.non_monophyletic_inc_sed_to_desIds[inc_sed_taxon->get_ott_id()] = inc_sed_taxon->get_data().desIds;
                break;
            }
        }
    }
    LOG(DEBUG) << " # non monophyletic incertae_sedis is " << unprune_stats.non_monophyletic_inc_sed_to_desIds.size();
}

// last step in unpruneTaxa
void unprunePreppedInputs(Tree_t & taxonomy,
                          Tree_t & solution,
                          const OttIdSet & incertae_sedis_ids,
                          const map<long, childParPair> & ott_to_tax,
                          UnpruneStats & unprune_stats) {
    map<long, Node_t*> ott_to_sol;
    const auto snVec = all_nodes(solution);
    // postorder walk over solution. Every time we find a node assigned to a taxon
    //  we augment the slice of the tree that is rooted at that node (and is the
    //  subtree that is cut at the deepest taxonomic node)
    // Note that the desIds field of soln is only filled in for slices. Once a 
    //    node with an OTT Id is hit, that ID will represent the entire clade for deeper
    //    nodes.
    auto & inc_sed_map = unprune_stats.inc_sed_taxon_to_sampled_tips;
    for (auto nd: snVec){
        if (nd->isTip()) {
            unprune_stats.startNumSolnLeaves += 1;
            if (!nd->has_ott_id()) {
                throw OTCError() << "Tip "<< nd->get_name() << " in solution lacks an OTT ID";
            }
            if (not ott_to_tax.count(nd->get_ott_id())) {
                throw OTCError() << "OttId "<< nd->get_ott_id() << " not in taxonomy!";
            }
        } else {
            if (nd->isOutDegreeOneNode()) {
                unprune_stats.startNumSolnMonotypicInternals += 1;
            } else {
                unprune_stats.startNumSolnForkingInternals += 1;
            }
        }
        auto p = nd->getParent();
        if (nd->has_ott_id()){
            auto ott_id = nd->get_ott_id();
            ott_to_sol[ott_id] = nd;
            auto & nd_di = nd->get_data().desIds;
            auto taxon_node = ott_to_tax.at(ott_id).first;
            if (nd->isTip()) {
                nd_di.insert(ott_id);
            }
            if (!taxon_node->isTip()) {
                unprune_stats.startNumNamedInternalsInSoln += 1;
                unpruneTaxaForSubtree(nd, ott_to_tax, ott_to_sol, unprune_stats);
                if (p) {
                    assert(p == nd->getParent());
                }
            } else {
                auto ismit = inc_sed_map.find(taxon_node);
                if (ismit != inc_sed_map.end() && !ismit->second.empty()) {
                    ismit->second.clear(); // flag this as a finished inc. sed. taxon because it is a tip
                }
                nd_di.insert(nd->get_ott_id());
            }
            if (p) {
                auto & pd = p->get_data();
                pd.desIds.insert(nd_di.begin(), nd_di.end());
            }
        } else {
            if (!nd->isTip()) {
                unprune_stats.numUnnamedInternalsInSoln += 1;
            }
            if (p) {
                const auto & d = nd->get_data().desIds;
                p->get_data().desIds.insert(d.begin(), d.end());
            }
        }
    }

    // Currently, we don't do the correct bookkeeping to get an accurate
    //    lost taxon map for non-monophyletic incertae sedis taxa during construction
    //    of the unpruned tree.
    //   So we'll just replace those entries here, so they are correct on exit.
    ott_to_sol.clear();
    for (auto n : iter_leaf(solution)) {
        ott_to_sol[n->get_ott_id()] = n;
    }
    const auto & id_to_fix_map = unprune_stats.non_monophyletic_inc_sed_to_desIds;
    for (auto pit : id_to_fix_map) {
        OttId broken_ott_id = pit.first;
        const auto & tip_ids = pit.second;
        fixLostTaxonMap(pit.first, pit.second, ott_to_sol, unprune_stats);
    }
}

// Note that this function will steal some nodes from `taxonomy`
//  so, on output that will no longer be a valid tree. This should
//  be fine if both taxonomy and solution go out of scope together.
template<typename T>
void unpruneTaxa(T & taxonomy,
                 T & solution,
                 UnpruneStats & unprune_stats) {
    const OttIdSet & incertae_sedis_ids = unprune_stats.incertae_sedis_ids;
    unprune_stats.out_degree_many1 = n_internal_out_degree_many(taxonomy);
    // 1. First, index taxonomy by OttId.
    map<long, childParPair> ott_to_tax;
    indexNodesByOttId(taxonomy, ott_to_tax, unprune_stats);
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
    // Note that at this point the desIds field of the taxonomy will be used for tracing
    //    parts of the solution. It will no longer be used as a flag for incertae sedis
    //    but the taxonomic nodes that were incertae sedis have been recorded in 
    //    unprune_stats.inc_sed_taxon_to_sampled_tips and the nonexcluded_ids sets
    for (auto nd: iter_post(taxonomy)) {
        nd->get_data().desIds.clear();
    }
    for (auto el : unprune_stats.inc_sed_taxon_to_sampled_tips) {
        LOG(DEBUG) << " unprune_stats.inc_sed_taxon_to_sampled_tips key " << el.first->get_ott_id() << " has " << el.second.size() << " sampled tips";
    }
    // Unprune, collecting stats
    unprunePreppedInputs(taxonomy, solution, incertae_sedis_ids, ott_to_tax, unprune_stats);
}

void reportStats(const UnpruneStats & unprune_stats, ostream * statsStreamPtr) {
    const OttIdSet & monotypicOttIds = unprune_stats.monotypic_ott_ids;
    const LostTaxonMap & ltm = unprune_stats.lost_taxa;
    auto numTaxaMonotypicInternals = monotypicOttIds.size();
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
        input["num_incertae_sedis_taxa"] = unprune_stats.incertae_sedis_ids.size();
        output["num_taxa_rejected"] = numTaxaRejected;
        output["num_taxonomy_splits_rejected"] = numForkingTaxaRejected;
        output["num_taxonomy_internals_merged"] = unprune_stats.numInternalsMergedOnBackbone;
        output["num_solution_leaves_expanded"] = unprune_stats.numLeavesExpanded;
        output["num_leaves_added"] = numLeavesAdded;
        document["input"] = input;
        document["output"] = output;
        *statsStreamPtr << document.dump(1) << std::endl;
    }
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
    UnpruneStats unprune_stats(incertae_sedis_ids);
    unpruneTaxa(taxonomy, solution, unprune_stats);
    reportStats(unprune_stats, statsStreamPtr);
    nameUnamedNodes(solution);
    writeTreeAsNewick(std::cout, solution);
    std::cout << std::endl;
    if (!lostTaxaJSONFilename.empty()) {
        writeLostTaxa(lostTaxaJSONStream,
                      unprune_stats.lost_taxa,
                      unprune_stats.synonym_map);
        lostTaxaJSONStream.close();
    }
    if (statsStreamPtr) {
        statsJSONStream.close();
    }
    return 0;
}
