#include <algorithm>
#include <unordered_set>
#include <iterator>
#include <sstream>
#include <stack>
#include <tuple>
#include <deque>
#include <utility>
#include "json.hpp"
#include "otc/lost_taxon_info.h"
#include "otc/otcli.h"
#include "otc/node_naming.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"

using namespace otc;

static std::string lostTaxaJSONFilename;
static std::string statsJSONFilename;
static std::string incertaeSedisFilename;
using std::forward_as_tuple;
using std::piecewise_construct;



struct RTNodePartialDesSet;
using Node_t = RootedTreeNode<RTNodePartialDesSet>;

typedef LostTaxonDetails<Node_t> LostTaxInfo;
typedef std::map<OttId, LostTaxInfo > LostTaxonMap;

struct RTNodePartialDesSet {
    std::set<long> des_ids;
    OttIdSet * nonexcluded_ids = nullptr; // union of IDs that descend from any child of an ancestor marked as incertae sedis
    OttId smallest_child  = 0;
    int tax_level = -1;
    // set for descendants of incertae sedis taxa to indicate which clade hold the taxon
    Node_t * repr_in_solution_by = nullptr; 
};

struct RTTreeIncertaeSedisHolder {
    std::list<OttIdSet> owned_id_sets; // used for memory management.
};

using Tree_t = RootedTree<RTNodePartialDesSet, RTTreeIncertaeSedisHolder>;


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



// adds ottId to the des_ids field of every node from firstNd down to ancAndLast (inclusive)
void addToDesIdsForAnc(long ottId, Node_t *firstNd, Node_t *ancAndLast, const map<long, childParPair> & ott2tax) {
    LOG(DEBUG) << "Adding " << ottId << " to des_ids for all taxa from " << firstNd->get_ott_id() << " to " << ancAndLast->get_ott_id();
    assert(firstNd);
    while (true) {
        firstNd->get_data().des_ids.insert(ottId);
        //LOG(DEBUG) << "addToDesIdsForAnc hit " << firstNd->get_ott_id();
        //LOG(DEBUG) << "adding " << ottId << " to node " << firstNd->get_name(); db_write_ott_id_set(" its des_ids", firstNd->get_data().des_ids);
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
    assert(!curr->get_data().des_ids.empty());
    assert(curr->get_data().repr_in_solution_by == nullptr);
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
                    if (!cd.des_ids.empty() && cd.repr_in_solution_by == nullptr) {
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
    solnChild->replace_this_node(nn);
    nn->set_name(taxon->get_name());
    nn->set_ott_id(taxon->get_ott_id());
    taxon->get_data().repr_in_solution_by = nn;
    nn->add_child(solnChild);
    nn->get_data().des_ids = solnChild->get_data().des_ids;
    return nn;
}

template <typename N>
N * addParentAndMoveUnsampledTaxChildren(N * taxon, N * solnChild) {
    auto * nn = introduceMonotypicParent(taxon, solnChild);
    moveUnsampledChildren(taxon, nn);
    return nn;
}

// Moves all of the children of taxon that have empty des_ids fields
//  to be children of solnNode.
template <typename N>
void moveUnsampledChildren(N * taxon, N * solnNode, std::set<N*> *v) {
    assert(taxon);
    assert(solnNode);
    const auto children = all_children(taxon);
    for (auto child : children) {
        auto & cd = child->get_data();
        if (cd.des_ids.empty() && cd.repr_in_solution_by == nullptr) {
            //LOG(DEBUG) << "moveUnsampledChildren moving " << child->get_ott_id();
            child->detach_this_node();
            solnNode->add_child(child);
            if (v) {
                v->insert(child);
            }
        }
    }
}


template <typename N>
N * special_bisect_with_new_child(N * taxNode, N * currSolnNd, bool moveUnsamp) {
    N * nt = bisect_branch_with_new_child(currSolnNd);
    nt->set_name(taxNode->get_name());
    nt->set_ott_id(taxNode->get_ott_id());
    nt->get_data().des_ids = currSolnNd->get_data().des_ids;
    if (moveUnsamp) {
        moveUnsampledChildren(taxNode, nt);
    }
    return nt;
}

// Looks through the children of `currSolnNd` If one has all the des_ids that are
//    flagged by taxon_node, it returns that node. Otherwise nullptr.
// Used for starting from an ancestor of the taxa in taxon_node and moving one
//    step closer to the MRCA of those taxa
Node_t * find_single_child_with_all_marked_taxa(Node_t * currSolnNd, const Node_t * taxon_node) {
    db_write_ott_id_set(" find_single_child_with_all_marked_taxa taxon des_ids", taxon_node->get_data().des_ids);
    db_write_ott_id_set(" find_single_child_with_all_marked_taxa currSolnNd des_ids", currSolnNd->get_data().des_ids);
    const auto & taxon_data = taxon_node->get_data();
    const auto & taxon_des = taxon_data.des_ids;
    Node_t * nextSolnNd = nullptr;
    for (auto c : iter_child(*currSolnNd)) {
        const auto & cdesId = c->get_data().des_ids;
        if (!are_disjoint(cdesId, taxon_des)) {
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
                                  const map<long, childParPair> & ott2tax,
                                  map<N *, set<N*> > & non_mono_to_register) {
    assert(higherTaxonNd);
    assert(rootSolnNd);
    LOG(DEBUG) << "higherTaxonNd->name = " << higherTaxonNd->get_name();
    const auto & taxData = higherTaxonNd->get_data();
    const auto & taxDes = taxData.des_ids;
    db_write_ott_id_set(" higherTaxonNd des_ids", taxDes);
    assert(taxDes.size() > 0);
    assert(taxData.nonexcluded_ids != nullptr);
    const auto & nonexcluded_ids = *taxData.nonexcluded_ids;
    N * currSolnNd = rootSolnNd;
    N * tipmostMRCA = nullptr;
    // the root of the solution, must have all of the constituent taxa
    assert(is_subset(taxDes, rootSolnNd->get_data().des_ids));
    // here we walk tipward to find the MRCA of the taxa in taxDes
    while (true) {
        N * nextSolnNd = find_single_child_with_all_marked_taxa(currSolnNd, higherTaxonNd);
        if (nextSolnNd == nullptr) {
            break;
        }
        currSolnNd = nextSolnNd;
    }
    if (unprune_stats.non_monophyletic_inc_sed.count(higherTaxonNd) > 0) {
        moveUnsampledChildren(higherTaxonNd, currSolnNd);
    }
    assert(currSolnNd != nullptr);
    const auto & solDes = currSolnNd->get_data().des_ids;
    assert(is_subset(taxDes, solDes));
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
        std::set<N *> v;
        moveUnsampledChildren(higherTaxonNd, currSolnNd, &v);
        if (recordInLTM) {
            for (auto c : v) {
                c->get_data().des_ids.clear();
                c->get_data().des_ids.insert(c->get_ott_id());
            }
            non_mono_to_register[higherTaxonNd] = v;
        }
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
    if (higherTaxonNd->is_outdegree_one_node()) {
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
            add_new_node = is_anc_des_pair(higherTaxonNd, potential_des);
            if ((!add_new_node) && rootSolnNd == currSolnNd) {
                add_new_node = is_anc_des_pair(potential_des, higherTaxonNd);
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
                currSolnNd->set_name(higherTaxonNd->get_name());
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
    effTipTaxonNd->get_data().des_ids.insert(effectiveTipOttId);
    auto effTipSolnNd = ott2soln.at(effectiveTipOttId);
    if (inc_sed_internals.find(effTipTaxonNd) != inc_sed_internals.end()) {
        curr_slice_inc_sed_map[effTipTaxonNd].insert(TaxSolnNdPair(effTipTaxonNd, effTipSolnNd));
    }
    while (anc &&  inc_sed_internals.find(anc) != inc_sed_internals.end()) {
        //LOG(DEBUG) << "registerCoveringOfIncSed added " << effectiveTipOttId << " to " << anc->get_ott_id();
        anc->get_data().des_ids.insert(effectiveTipOttId);
        curr_slice_inc_sed_map[anc].insert(TaxSolnNdPair(effTipTaxonNd, effTipSolnNd));
        anc = ott2tax.at(anc->get_ott_id()).second;
    }
}

// Moves unsampled taxa from the taxanomy (by finding them in `ott2tax`) to a slice of the
//  solution tree that is rooted at `rootSolnNd`.
// Uses the `rootSolnNd->get_data().des_ids` to denote the set of leaves of this slice of the
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
    assert(!rootTaxonNd->is_tip());
    const OttIdSet * rootNonExclude = rootTaxonNd->get_data().nonexcluded_ids;
    assert(rootNonExclude != nullptr);
    // this will have the IDs for all of the include taxa for this slice of the tree
    // note that some may be "higher" taxa.
    auto & solnDesIds = rootSolnNd->get_data().des_ids;
    assert(!solnDesIds.empty());
    // des_ids fields of the solution tree are filled in by the caller before this function.
    //  here we fill in those fields for the taxonomy nodes.
    //Here we add the sampled IDs to desId fields ofthe relevant taxa. This is potentially
    //  confusing because get_data().des_ids will hold only the IDs of taxa that are present in the soln
    //  tree (even when we are talking about the taxonomy node's des_ids)
    // While we are walking through the "leaf" Ids for this tree slice, we'll also
    //  collect the leaf sets
    const auto & inc_sed_map = unprune_stats.inc_sed_taxon_to_sampled_tips;
    const auto & inc_sed_internals = unprune_stats.inc_sed_internals;
    set<N *> solnLeaves;
    set<N *> taxaLeaves;
    std::map<Node_t *, std::set<TaxSolnNdPair> > curr_slice_inc_sed_map;
    set<Node_t *> inc_sed_that_are_proper_children;
    db_write_ott_id_set(" effective tips: ", solnDesIds);
    db_write_ott_id_set(" rootNonExclude: ", *rootNonExclude);
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
        if (!is_inc_sed && effTipTaxonNd->is_tip()) {
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
                if (effTipTaxonNd->is_tip()) {
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
            if (effTipTaxonNd->get_data().repr_in_solution_by != nullptr) {
                assert(effTipTaxonNd->get_data().repr_in_solution_by->has_ott_id());
                queuedForFlagging.push_back(carriedIDRealId(effectiveTipOttId, effTipTaxonNd->get_data().repr_in_solution_by->get_ott_id()));
                LOG(DEBUG) << "checking effective tip " << effectiveTipOttId << " threaded through " << effTipTaxonNd->get_data().repr_in_solution_by->get_ott_id();
            } else {
                queuedForFlagging.push_back(carriedIDRealId(effectiveTipOttId, effectiveTipOttId));
            }
        }    
    }
    bool root_taxon_is_inc_sed = inc_sed_map.find(rootTaxonNd) != inc_sed_map.end();
    if (root_taxon_is_inc_sed) {
        inc_sed_that_are_proper_children.insert(rootTaxonNd);
    }
    // now we flag all of the des_ids for the proper
    for (auto cr_pair : queuedForFlagging) {
        auto realTipOttId = cr_pair.second;
        auto effTipTaxonNd = ott2tax.at(realTipOttId).first;
        effTipTaxonNd->get_data().des_ids.clear();
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
    db_write_ott_id_set("Incertae sedis taxa for the current slice...", db_inc_sed_ott_ids);
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
        if (l->is_tip()) {
            assert(l->has_ott_id());
            auto leafOttId = l->get_ott_id();
            auto taxonForLeaf = ott2tax.at(leafOttId).first;
            if (!taxonForLeaf->is_tip()) {
                numExpanded += 1;
                auto cvec = all_children(taxonForLeaf);
                for (auto c : cvec) {
                    if (c->get_data().repr_in_solution_by == nullptr) {
                        taxonForLeaf->remove_child(c);
                        c->get_data().repr_in_solution_by = l;
                        l->add_child(c);
                    }
                }
            }
        }
    }
    size_t numMerged = numExpanded;
    // Walk through the taxonomy in a preorder fashion, but only accumulate
    //  a vector of those nodes with at least 1 member of the `des_ids` field.
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
    map<N *, set<N *> > non_mono_to_register;
    for (auto higherTaxonNd : postOrderInTaxNd) {
        LOG(DEBUG) << "dealing with " << higherTaxonNd->get_ott_id();
        if (higherTaxonNd == rootTaxonNd) {
            moveUnsampledChildren(higherTaxonNd, rootSolnNd);
            numMerged += 1 ;
        } else {
            bool recordInLTM = inc_sed_that_are_proper_children.find(higherTaxonNd) == inc_sed_that_are_proper_children.end();
            numMerged += incorporateHigherTaxonNode(higherTaxonNd, rootSolnNd, nodesAddedForTaxa, unprune_stats, recordInLTM, ott2tax, non_mono_to_register);
        }
        LOG(DEBUG) << "repr_in_solution_by marked for " << higherTaxonNd->get_ott_id();
        higherTaxonNd->get_data().repr_in_solution_by = rootSolnNd;
    }
    for (auto is_it = inc_sed_to_deal_with_by_level.rbegin(); is_it != inc_sed_to_deal_with_by_level.rend(); ++is_it) {
        for (auto is_taxon_ptr : is_it->second) {
            LOG(DEBUG) << "dealing with Incertae sedis " << is_taxon_ptr->get_ott_id() << " from level " << is_it->first;
            if (is_taxon_ptr == rootTaxonNd) {
                moveUnsampledChildren(is_taxon_ptr, rootSolnNd);
                numMerged += 1 ;
            } else {
                numMerged += incorporateHigherTaxonNode(is_taxon_ptr, rootSolnNd, nodesAddedForTaxa, unprune_stats, false, ott2tax, non_mono_to_register);
            }
            LOG(DEBUG) << "repr_in_solution_by marked for " << is_taxon_ptr->get_ott_id();
            is_taxon_ptr->get_data().repr_in_solution_by = rootSolnNd;
        }
    }
    for (auto nmel : non_mono_to_register) {
        auto nm = nmel.first;
        const auto & moved_nodes = nmel.second;
        for (auto mtn : moved_nodes) {
            ott2soln[mtn->get_ott_id()] = mtn;
        }
        const auto moved_ids = nodes_to_ott_id_set(moved_nodes);
        auto & taxon_des = nm->get_data().des_ids;
        const auto * ni = nm->get_data().nonexcluded_ids;
        assert(ni != nullptr);
        nm->get_data().des_ids.insert(moved_ids.begin(), moved_ids.end());
        assert(nm->has_ott_id());
        auto nm_ott_id = nm->get_ott_id();
        unprune_stats.lost_taxa.emplace(piecewise_construct,
                                        forward_as_tuple(nm_ott_id),
                                        forward_as_tuple(nm_ott_id, taxon_des, *ni, ott2soln));
    }
    // for the sake of a low memory footprint, here we 
    //  clear out the des_ids fields for this slice of the tree.
    for (auto n : taxaLeaves) {
        n->get_data().repr_in_solution_by = rootSolnNd;
        n->get_data().des_ids.clear();
    }
    for (auto n : postOrderInTaxNd) {
        n->get_data().des_ids.clear();
    }
    set<N *> visited;
    for (auto n : solnLeaves) {
        auto c = n;
        while (visited.count(c) == 0 && c != rootSolnNd) {
            c->get_data().des_ids.clear();
            if (c != n) {
                visited.insert(c);
            }
            c = c->get_parent();
        }
    }
    for (auto scism_it: curr_slice_inc_sed_map) {
        auto & inc_sed_taxon = scism_it.first;
        inc_sed_taxon->get_data().des_ids.clear();
    }
    rootSolnNd->get_data().des_ids.clear();
    rootSolnNd->get_data().des_ids.insert(ottId);
    // here we need to update the altered inc_sed taxon maps...
    // All the ones finished in this slice can be deleted from maps...
    for (auto is_taxon_ptr : inc_sed_that_are_proper_children) {
        LOG(DEBUG) << "removing des " << is_taxon_ptr->get_ott_id() << " from inc_sed_taxon_to_sampled_tips";
        unprune_stats.inc_sed_taxon_to_sampled_tips.at(is_taxon_ptr).clear();
        if (is_taxon_ptr->get_parent()) {
            is_taxon_ptr->detach_this_node();
        }
    }
    for (auto is_it = inc_sed_to_deal_with_by_level.rbegin(); is_it != inc_sed_to_deal_with_by_level.rend(); ++is_it) {
        for (auto is_taxon_ptr : is_it->second) {
            LOG(DEBUG) << "removing " << is_taxon_ptr->get_ott_id() << " from inc_sed_taxon_to_sampled_tips";
            unprune_stats.inc_sed_taxon_to_sampled_tips.at(is_taxon_ptr).clear();
            if (is_taxon_ptr->get_parent()) {
                is_taxon_ptr->detach_this_node();
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
            rootSolnNd->get_data().des_ids.insert(spike_nd->get_ott_id());
            //flagIncSedAs
            while (true) {
                if (inc_sed_mapping_deeper.count(spike_nd) > 0) {
                    break;
                }
                LOG(DEBUG) << "repr_in_solution_by marked for " << spike_nd->get_ott_id();
                spike_nd->get_data().repr_in_solution_by = rootTaxonNd;
                spike_nd = spike_nd->get_parent();
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
    auto & owningList = taxonomy.get_data().owned_id_sets;
    auto taxonomy_root_ptr = taxonomy.get_root();
    owningList.push_back(OttIdSet());
    taxonomy_root_ptr->get_data().nonexcluded_ids = &(*owningList.rbegin());
    taxonomy_root_ptr->get_data().tax_level = 0;
    for (auto nd: iter_pre(taxonomy)) {
        if (nd->is_tip()) {
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
                    const auto & isc_di = isc->get_data().des_ids;
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
                    auto isc_di = isc->get_data().des_ids;
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
    // As a part of implementing the incertae sedis logic, we also fill the des_ids of the taxonomy.
    // This info is used to create the "non-excluded" sets, but it is expensive to do on the whole 
    //    taxonomy, so we just do it for subtrees rooted at an incertae sedis taxon...
    for (auto ist_id: incertae_sedis_ids) {
        //LOG(DEBUG) << ""
        if (ott_to_tax.count(ist_id) == 0) {
            throw OTCError() << "Incertae sedis ID (" << ist_id << ") not in taxonomy.";
        }
        LOG(DEBUG) << "Filling des_ids for incertae sedis " << ist_id;
        auto taxon = ott_to_tax.at(ist_id).first;
        for (auto nd: iter_post_n(*taxon)) {
            auto & di_ref = nd->get_data().des_ids;
            // empty check is just an optimization to avoid recalculating the results for a node.
            if (di_ref.empty()) {
                di_ref.insert(nd->get_ott_id());
                unprune_stats.inc_sed_internals.insert(nd);
                for (auto c : iter_child(*nd)) {
                    const auto & cdi_ref = c->get_data().des_ids;
                    di_ref.insert(cdi_ref.begin(), cdi_ref.end());
                }
            }
        }
        db_write_ott_id_set("the des_ids set is ", taxon->get_data().des_ids);
    }
}


void indexNodesByOttId(Tree_t & taxonomy,
                       map<long, childParPair> & ott_to_tax,
                       UnpruneStats & unprune_stats) {
    OttIdSet & moi = unprune_stats.monotypic_ott_ids;
    for (auto nd: iter_post(taxonomy)){
        if (nd->has_ott_id()) {
            const OttId taxon_ott_id = nd->get_ott_id();
            ott_to_tax[taxon_ott_id] = childParPair(nd, nd->get_parent());
            if (nd->is_outdegree_one_node()) {
                moi.insert(nd->get_ott_id());
            } else if (nd->is_tip()) {
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
            db_write_ott_id_set("target des_ids", taxon_node->get_data().des_ids);
        }
        if (taxon_node->get_data().des_ids.empty()) {
            continue;    
        }
        unprune_stats.inc_sed_soln_tips.insert(snd);
        auto anc = taxon_node->get_parent();
        to_tips_map[taxon_node].insert(snd);
        while (anc && !(anc->get_data().des_ids.empty())) {
            to_tips_map[anc].insert(snd);
            anc = anc->get_parent();
        }
    }
    // db
    for (auto x : to_tips_map) {
        LOG(DEBUG) << "findIncertaeSedisInSolution recorded entry in inc_sed_taxon_to_sampled_tips for " << x.first->get_ott_id();    
    }
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
        const auto taxon_des = inc_sed_taxon->get_data().des_ids;
        assert(nonexc != nullptr);
        auto soln_mrca = find_mrca_via_traversing(sampled);
        for (auto t : iter_leaf_n(*soln_mrca)) {
            LOG(DEBUG) << " visiting leaf " << t->get_ott_id() << " to see if it breaks " << inc_sed_taxon->get_ott_id();
            if (sampled.count(t) == 0 && 0 == nonexc->count(t->get_ott_id())) {
                LOG(DEBUG) << "Inc. sed. " << inc_sed_taxon->get_ott_id() << " is not monophyletic.";
                unprune_stats.non_monophyletic_inc_sed[inc_sed_taxon] = soln_mrca;
                unprune_stats.non_monophyletic_inc_sed_to_desIds[inc_sed_taxon->get_ott_id()] = inc_sed_taxon->get_data().des_ids;
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
    // Note that the des_ids field of soln is only filled in for slices. Once a 
    //    node with an OTT Id is hit, that ID will represent the entire clade for deeper
    //    nodes.
    auto & inc_sed_map = unprune_stats.inc_sed_taxon_to_sampled_tips;
    for (auto nd: snVec){
        if (nd->is_tip()) {
            unprune_stats.startNumSolnLeaves += 1;
            if (!nd->has_ott_id()) {
                throw OTCError() << "Tip "<< nd->get_name() << " in solution lacks an OTT ID";
            }
            if (not ott_to_tax.count(nd->get_ott_id())) {
                throw OTCError() << "OttId "<< nd->get_ott_id() << " not in taxonomy!";
            }
        } else {
            if (nd->is_outdegree_one_node()) {
                unprune_stats.startNumSolnMonotypicInternals += 1;
            } else {
                unprune_stats.startNumSolnForkingInternals += 1;
            }
        }
        auto p = nd->get_parent();
        if (nd->has_ott_id()){
            auto ott_id = nd->get_ott_id();
            ott_to_sol[ott_id] = nd;
            auto & nd_di = nd->get_data().des_ids;
            auto taxon_node = ott_to_tax.at(ott_id).first;
            if (nd->is_tip()) {
                nd_di.insert(ott_id);
            }
            if (!taxon_node->is_tip()) {
                unprune_stats.startNumNamedInternalsInSoln += 1;
                unpruneTaxaForSubtree(nd, ott_to_tax, ott_to_sol, unprune_stats);
                if (p) {
                    assert(p == nd->get_parent());
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
                pd.des_ids.insert(nd_di.begin(), nd_di.end());
            }
        } else {
            if (!nd->is_tip()) {
                unprune_stats.numUnnamedInternalsInSoln += 1;
            }
            if (p) {
                const auto & d = nd->get_data().des_ids;
                p->get_data().des_ids.insert(d.begin(), d.end());
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
    auto & lost_taxa = unprune_stats.lost_taxa;
    for (auto pit : id_to_fix_map) {
        OttId broken_ott_id = pit.first;
        auto taxon_node = ott_to_tax.at(broken_ott_id).first;
        auto nonexcluded_ids = taxon_node->get_data().nonexcluded_ids;
        assert(nonexcluded_ids != nullptr);
        lost_taxa.emplace(piecewise_construct,
                          forward_as_tuple(broken_ott_id),
                          forward_as_tuple(broken_ott_id,
                                           pit.second,
                                           *nonexcluded_ids,
                                           ott_to_sol));
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
    // filling the nonexcluded IDs (next step) needs the des_ids filled for the parts of 
    //    the taxonomy that descend from an incertae sedis taxon.
    // Note that we'll use empty des_ids of a taxonomy node as a flag that that taxon is 
    //    not a descendant of an incertae sedis taxon.
    fillDesIdsForIncertaeSedisOnly(ott_to_tax, incertae_sedis_ids, unprune_stats);
    // Now we fill in the unique non-excluded Ids. A non-excluded ID for an internal
    //    node y is the ID of a taxon that is a descendant of some incertae sedis taxon x
    //    where x is child of one of the ancestors of the node y.
    fillNonexcludeIDFields(taxonomy, incertae_sedis_ids);
    // record tips in the sample that are incertae sedis descendants
    findIncertaeSedisInSolution(taxonomy, solution, ott_to_tax, unprune_stats);
    // note which ones can't be in the solution
    findBrokenIncertaeSedisDescendants(taxonomy, solution, ott_to_tax, unprune_stats);
    // Note that at this point the des_ids field of the taxonomy will be used for tracing
    //    parts of the solution. It will no longer be used as a flag for incertae sedis
    //    but the taxonomic nodes that were incertae sedis have been recorded in 
    //    unprune_stats.inc_sed_taxon_to_sampled_tips and the nonexcluded_ids sets
    for (auto nd: iter_post(taxonomy)) {
        nd->get_data().des_ids.clear();
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
    otCLI.add_flag('j',
                  "Produce a JSON file that notes where the descendants of unincluded appear on the tree",
                  handleLostTaxaJSON,
                  true);
    otCLI.add_flag('s',
                  "Produce a JSON file with stats about the inputs and outputs",
                  handleStatsJSON,
                  true);
    otCLI.add_flag('i',
                  "Optional list of IDs of tree in the exemplified taxonomy that are incertae sedis",
                  handleIncertaeSedis,
                  true);
    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};
    if (argc < 2) {
        throw OTCError("No trees provided!");
    }
    if (tree_processing_main<Tree_t>(otCLI, argc, argv, get, nullptr, 1)) {
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
    name_unnamed_nodes(solution);
    write_tree_as_newick(std::cout, solution);
    std::cout << std::endl;
    if (!lostTaxaJSONFilename.empty()) {
        const LostTaxonMap & ltm = unprune_stats.lost_taxa;
        write_lost_taxa(lostTaxaJSONStream, ltm, unprune_stats.synonym_map);
        lostTaxaJSONStream.close();
    }
    if (statsStreamPtr) {
        statsJSONStream.close();
    }
    return 0;
}
