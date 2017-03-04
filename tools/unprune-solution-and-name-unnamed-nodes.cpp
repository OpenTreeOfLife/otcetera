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
using std::forward_as_tuple;
using std::piecewise_construct;

static std::string lostTaxaJSONFilename;
static std::string statsJSONFilename;
static std::string incertaeSedisFilename;

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
void move_unsampled_tax_children(N * taxon, N * solution_node, std::set<N*> *v = nullptr);

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



// adds ott_id to the des_ids field of every node from first_node down to anc_and_last_node (inclusive)
void add_to_des_ids_for_anc(long ott_id,
                            Node_t *first_node,
                            Node_t *anc_and_last_node,
                            const map<long, childParPair> & ott_to_tax) {
    assert(first_node);
    while (true) {
        first_node->get_data().des_ids.insert(ott_id);
        if (first_node == anc_and_last_node) {
            return;
        }
        auto fid = first_node->get_ott_id();
        first_node = ott_to_tax.at(fid).second; // parent in unaltered taxonomy.
        // If this triggers, then anc_and_last_node is not an ancestor of first_node
        assert(first_node);
    }
}

// Returns a vector of all of the nodes in the subtree rooted at `root`
// that are:
//    1. rootward of the boundary nodes,
//    2. has a non-empty des_ids field, and 
//    3. has a nullptr for the repr_in_solution_by field
// The nodes are returned in a preorder ordering.
//
// Prior to calling this function, the des_ids nodes for a slice of the tree are
//    filled in from the included tips down to the `root`. Nodes that have empty
//    des_ids are nodes that are not sampled in the solution.  Nodes that have
//    a non-NULL repr_in_solution_by field were processed in a previous slice
//    and are acting as tips in the current slice - so they act like boundary nodes.
// The taxa that are tips of the slice are passed in as the boundaries. 
// Thus the returned vector is a preorder traversal of the induced tree for the slice
//    that does not include the tips of that tree.  Thus, it is a vector of higher
//    taxa that could be added to this slice of the tree.
template <typename N>
vector<N *> preorder_below_boundaries(N * root,
                                      const set<N *> & boundaries) {
    set<N *> seen;
    vector<N *> r;
    N * curr = root;
    assert(!curr->get_data().des_ids.empty());
    assert(curr->get_data().repr_in_solution_by == nullptr);
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

// introduces and returns a node in the solution as the parent of `solution_child`
// The name and ott_id are obtained from `taxon` and `taxon` has its
//    repr_in_solution_by field set to the returned node to indicate that 
//    the taxon has been added to the solution. 
template <typename N>
N * introduce_monotypic_parent(N * taxon, N * solution_child) {
    assert(taxon);
    assert(solution_child);
    N * nn = new N(nullptr);
    solution_child->replace_this_node(nn);
    nn->set_name(taxon->get_name());
    nn->set_ott_id(taxon->get_ott_id());
    taxon->get_data().repr_in_solution_by = nn;
    nn->add_child(solution_child);
    nn->get_data().des_ids = solution_child->get_data().des_ids;
    return nn;
}


// Given a node in the tree that represents the mrca of a `taxon`, but which
//    has already been assigned a distinct OTT Id, this function is used to
//    incorporate a node with the ID of `taxon` and all of its unsampled children
//    as the parent to `solution_child`
template <typename N>
N * add_parent_and_move_unsampled_tax_children(N * taxon, N * solution_child) {
    auto * nn = introduce_monotypic_parent(taxon, solution_child);
    move_unsampled_tax_children(taxon, nn);
    return nn;
}

// Moves all of the children of `taxon` that have empty des_ids fields
//  to be children of `solution_node`.
// If `v` is non-null, then the moved nodes are inserted into the set that 
//    `v` points to.
// Note that this function moves nodes from the taxonomy tree into the solution tree.
template <typename N>
void move_unsampled_tax_children(N * taxon, N * solution_node, std::set<N*> *v) {
    assert(taxon);
    assert(solution_node);
    const auto children = all_children(taxon);
    for (auto child : children) {
        auto & cd = child->get_data();
        if (cd.des_ids.empty() && cd.repr_in_solution_by == nullptr) {
            child->detach_this_node();
            solution_node->add_child(child);
            if (v) {
                v->insert(child);
            }
        }
    }
}

// This function is used when we want to add a new node for `taxon`, but it is introducing
//     a node at the root of the  current slice (`curr_solution_node`).
// Because we have aliases to that node, we need to add the taxon as a child of 
//    `curr_solution_node`
// This function renders `curr_solution_node` by adding a new node as a child and then
//    transferring the original children to that node. It gets its ID and des_ids
//    from `taxon`. If `move_unsamp` is true, then it also gets the unsampled taxa
//    that are children of `taxon via move_unsampled_tax_children.
template <typename N>
N * special_bisect_with_new_child(N * taxon, N * curr_solution_node, bool move_unsamp) {
    N * nt = bisect_branch_with_new_child(curr_solution_node);
    nt->set_name(taxon->get_name());
    nt->set_ott_id(taxon->get_ott_id());
    nt->get_data().des_ids = curr_solution_node->get_data().des_ids;
    if (move_unsamp) {
        move_unsampled_tax_children(taxon, nt);
    }
    return nt;
}

// Looks through the children of `curr_solution_node` If one has all the des_ids that are
//    flagged by taxon_node, it returns that node. Otherwise nullptr.
// Used for starting from an ancestor of the taxa in taxon_node and moving one
//    step closer to the MRCA of those taxa
Node_t * find_single_child_with_all_marked_taxa(Node_t * curr_solution_node,
                                                const Node_t * taxon_node) {
    const auto & taxon_data = taxon_node->get_data();
    const auto & taxon_des = taxon_data.des_ids;
    Node_t * next_solution_node = nullptr;
    for (auto c : iter_child(*curr_solution_node)) {
        const auto & cdesId = c->get_data().des_ids;
        if (!are_disjoint(cdesId, taxon_des)) {
            if (next_solution_node == nullptr) {
                next_solution_node = c;
            } else {
                return nullptr;
            }
        }
    }
    return next_solution_node;
}

// returns the number of nodes on the backbone that are merged with this higher taxon
//      because the taxon was compatible (the same as) the node - will be either 0 or 1.
template <typename N>
size_t incorporate_higher_taxon(N* taxon,
                                N* root_solution_node,
                                set<N *> & nodes_add_for_taxa,
                                UnpruneStats & unprune_stats,
                                bool record_in_ltm,
                                const map<long, childParPair> & ott_to_tax,
                                map<N *, set<N*> > & non_mono_to_register) {
    assert(taxon);
    assert(root_solution_node);
    const auto & taxData = taxon->get_data();
    const auto & taxDes = taxData.des_ids;
    assert(taxDes.size() > 0);
    assert(taxData.nonexcluded_ids != nullptr);
    const auto & nonexcluded_ids = *taxData.nonexcluded_ids;
    N * curr_solution_node = root_solution_node;
    N * tipmostMRCA = nullptr;
    // the root of the solution, must have all of the constituent taxa
    assert(is_subset(taxDes, root_solution_node->get_data().des_ids));
    // here we walk tipward to find the MRCA of the taxa in taxDes
    while (true) {
        N * next_solution_node = find_single_child_with_all_marked_taxa(curr_solution_node, taxon);
        if (next_solution_node == nullptr) {
            break;
        }
        curr_solution_node = next_solution_node;
    }
    if (unprune_stats.non_monophyletic_inc_sed.count(taxon) > 0) {
        move_unsampled_tax_children(taxon, curr_solution_node);
    }
    assert(curr_solution_node != nullptr);
    const auto & solDes = curr_solution_node->get_data().des_ids;
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
        std::set<N *> v;
        move_unsampled_tax_children(taxon, curr_solution_node, &v);
        if (record_in_ltm) {
            for (auto c : v) {
                c->get_data().des_ids.clear();
                c->get_data().des_ids.insert(c->get_ott_id());
            }
            non_mono_to_register[taxon] = v;
        }
        return 0U;
    }
    // This solution node has the same descendant set of this taxon.
    // It could either be the only child of the taxon (if the taxon) is 
    //    monotypic, or it could be labelled with this OTT ID.
    N * nt = nullptr;
    // the easiest case is if the solution node has the same taxon ID
    //    just move the unsampled children.
    if (curr_solution_node->has_ott_id() && curr_solution_node->get_ott_id() == taxon->get_ott_id()) {
        move_unsampled_tax_children(taxon, curr_solution_node);
        return 1U;
    }
    // The monotypic case is also easy... add a monotypic node to the solution.
    if (taxon->is_outdegree_one_node()) {
        // careful to keep the root of the subproblem the same node structure, which
        //    may require swapping data to make the root the new deeper monotypic taxon....
        if (curr_solution_node == root_solution_node) {
            nt = special_bisect_with_new_child(taxon, curr_solution_node, false);
        } else {
            nt = introduce_monotypic_parent(taxon, curr_solution_node);
        }
    } else {
        // In the world of incertae sedis, we need to worry about the case of the node
        //    already having a label, but not being a descendant of higher taxon.
        // We only want to add a new node to the solution if the solution node that
        //    we have found is associated with a descendant taxon of
        bool add_new_node = false;
        if (curr_solution_node->has_ott_id()) {
            // we have to check the taxonomy tree to find out...
            auto potential_des = ott_to_tax.at(curr_solution_node->get_ott_id()).first;
            add_new_node = is_anc_des_pair(taxon, potential_des);
            if ((!add_new_node) && root_solution_node == curr_solution_node) {
                add_new_node = is_anc_des_pair(potential_des, taxon);
            }
        }
        if (add_new_node) {
            if (curr_solution_node == root_solution_node) {
                nt = special_bisect_with_new_child(taxon, curr_solution_node, true);
            } else {
                nt = add_parent_and_move_unsampled_tax_children(taxon, curr_solution_node);
            }
        } else {
            move_unsampled_tax_children(taxon, curr_solution_node);
            if (curr_solution_node->has_ott_id()) {
                unprune_stats.synonym_map[curr_solution_node->get_ott_id()].insert(taxon->get_ott_id());
            } else {
                curr_solution_node->set_name(taxon->get_name());
                curr_solution_node->set_ott_id(taxon->get_ott_id());
            }
        }
    }
    if (nt != nullptr) {
        nodes_add_for_taxa.insert(nt);
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
                              const map<long, childParPair> & ott_to_tax) {
    const auto & inc_sed_internals = unprune_stats.inc_sed_internals;
    auto anc = ott_to_tax.at(effectiveTipOttId).second;
    effTipTaxonNd->get_data().des_ids.insert(effectiveTipOttId);
    auto effTipSolnNd = ott2soln.at(effectiveTipOttId);
    if (inc_sed_internals.find(effTipTaxonNd) != inc_sed_internals.end()) {
        curr_slice_inc_sed_map[effTipTaxonNd].insert(TaxSolnNdPair(effTipTaxonNd, effTipSolnNd));
    }
    while (anc &&  inc_sed_internals.find(anc) != inc_sed_internals.end()) {
        anc->get_data().des_ids.insert(effectiveTipOttId);
        curr_slice_inc_sed_map[anc].insert(TaxSolnNdPair(effTipTaxonNd, effTipSolnNd));
        anc = ott_to_tax.at(anc->get_ott_id()).second;
    }
}

// Moves unsampled taxa from the taxanomy (by finding them in `ott_to_tax`) to a slice of the
//  solution tree that is rooted at `root_solution_node`.
// Uses the `root_solution_node->get_data().des_ids` to denote the set of leaves of this slice of the
//  tree. On exit this will be the ID of this clade and the IDs of any uncoalesced incertae sedis 
//    taxa that are in this slice.
// If the taxon conflicts with the solution, then the unsampled taxa are attached at the
//  MRCA of the taxon in the solution, and an entry is added mapping the taxon's OTT Id
//  to a LostTaxonDetails object that explains where the elements of the taxon are located.
// Returns the set of IDs that should  
template <typename N>
void unpruneTaxaForSubtree(N *root_solution_node,
                           const map<long, childParPair> & ott_to_tax, 
                           map<long, N*> & ott2soln, 
                           UnpruneStats & unprune_stats) {
    assert(root_solution_node);
    assert(root_solution_node->has_ott_id());
    const auto ott_id = root_solution_node->get_ott_id();
    assert(ott2soln.at(ott_id) == root_solution_node);
    N * rootTaxonNd = ott_to_tax.at(ott_id).first;
    assert(rootTaxonNd != nullptr);
    assert(!rootTaxonNd->is_tip());
    const OttIdSet * rootNonExclude = rootTaxonNd->get_data().nonexcluded_ids;
    assert(rootNonExclude != nullptr);
    // this will have the IDs for all of the include taxa for this slice of the tree
    // note that some may be "higher" taxa.
    auto & solnDesIds = root_solution_node->get_data().des_ids;
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
    list<Node_t *> effTipTaxa;
    typedef pair<OttId, OttId> CarriedIDRealId;
    list<CarriedIDRealId > queuedForFlagging;
    for (auto effectiveTipOttId : solnDesIds) {
        auto effTipTaxonNd = ott_to_tax.at(effectiveTipOttId).first;
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
                // this "tip" taxon is nonexcluded incertae sedis taxon from deeper in tree...
                if (!is_tip_inc_sed && is_map_it->second.size() > 0) {
                    throw OTCError() << "Incertae sedis taxon " << effectiveTipOttId 
                                     << " is a tip for slice " << ott_id 
                                     << ", but is still in the inc_sed_taxon_to_sampled_tips map!";
                }
                registerCoveringOfIncSed(effectiveTipOttId, 
                                         effTipTaxonNd,
                                         ott2soln, 
                                         unprune_stats,
                                         curr_slice_inc_sed_map,
                                         ott_to_tax);
                if (effTipTaxonNd->is_tip()) {
                    taxaLeaves.insert(effTipTaxonNd);
                }
            } else {
                is_proper_des = true;
            }
        } else {
            is_proper_des = true;
        }
        if (is_proper_des) {
            if (is_inc_sed && !is_tip_inc_sed) {
                inc_sed_that_are_proper_children.insert(effTipTaxonNd);
            }
            if (effTipTaxonNd->get_data().repr_in_solution_by != nullptr) {
                assert(effTipTaxonNd->get_data().repr_in_solution_by->has_ott_id());
                auto cff = CarriedIDRealId(effectiveTipOttId,
                                           effTipTaxonNd->get_data().repr_in_solution_by->get_ott_id());
                queuedForFlagging.push_back(cff);
            } else {
                queuedForFlagging.push_back(CarriedIDRealId(effectiveTipOttId, effectiveTipOttId));
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
        auto effTipTaxonNd = ott_to_tax.at(realTipOttId).first;
        effTipTaxonNd->get_data().des_ids.clear();
    }
    for (auto cr_pair : queuedForFlagging) {
        auto carriedTipId = cr_pair.first;
        auto threadedThroughOttId = cr_pair.second;
        Node_t * effTipTaxonNd = ott_to_tax.at(carriedTipId).first;
        add_to_des_ids_for_anc(carriedTipId, effTipTaxonNd, rootTaxonNd, ott_to_tax);
        solnLeaves.insert(ott2soln.at(threadedThroughOttId));
        taxaLeaves.insert(effTipTaxonNd);
    }   
    // need to process inc. sed. splits in reverse order by level, to mimic the postorder sweep
    map<int, list<Node_t *> > inc_sed_to_deal_with_by_level;
    set<Node_t *> inc_sed_mapping_deeper;

    // Now we can use curr_slice_inc_sed_map to figure out which inc. sedis. taxa will be dealt with in this slice;
    for (auto scism_it: curr_slice_inc_sed_map) {
        auto & inc_sed_taxon = scism_it.first;
        auto & included_leaf_pair_set = scism_it.second;
        auto gsit = inc_sed_map.find(inc_sed_taxon);
        if (gsit == inc_sed_map.end()) {
            throw OTCError() << "Incertae sedis taxon " << inc_sed_taxon->get_ott_id() << " included in slice " << ott_id << ", but was not found in the inc_sed_taxon_to_sampled_tips map!";
        }
        auto & full_soln_set = gsit->second;
        if (full_soln_set.size() == 0) {
            continue; // we don't need to worry about incertae sedis taxa with no descendants stored. these have been dealt with
        }
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
            auto taxonForLeaf = ott_to_tax.at(leafOttId).first;
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
    auto postOrderInTaxNd = preorder_below_boundaries(rootTaxonNd, taxaLeaves);
    std::reverse(postOrderInTaxNd.begin(), postOrderInTaxNd.end());
    // Here we will add the relevant higher (non-leaf) taxa from the taxonomy to the 
    //  solution. It is crucial that we do this in postorder because if we have a series
    //  of nodes to introduce along a branch, the incorporate_higher_taxon function
    //  will add them below the attachment node (so postorder will assure that they 
    //  are correctly added in the tip->root orientation).
    set<N *> nodes_add_for_taxa;
    map<N *, set<N *> > non_mono_to_register;
    for (auto taxon : postOrderInTaxNd) {
        if (taxon == rootTaxonNd) {
            move_unsampled_tax_children(taxon, root_solution_node);
            numMerged += 1 ;
        } else {
            bool record_in_ltm = inc_sed_that_are_proper_children.find(taxon) == inc_sed_that_are_proper_children.end();
            numMerged += incorporate_higher_taxon(taxon, root_solution_node, nodes_add_for_taxa, unprune_stats, record_in_ltm, ott_to_tax, non_mono_to_register);
        }
        taxon->get_data().repr_in_solution_by = root_solution_node;
    }
    for (auto is_it = inc_sed_to_deal_with_by_level.rbegin(); is_it != inc_sed_to_deal_with_by_level.rend(); ++is_it) {
        for (auto is_taxon_ptr : is_it->second) {
            if (is_taxon_ptr == rootTaxonNd) {
                move_unsampled_tax_children(is_taxon_ptr, root_solution_node);
                numMerged += 1 ;
            } else {
                numMerged += incorporate_higher_taxon(is_taxon_ptr, root_solution_node, nodes_add_for_taxa, unprune_stats, false, ott_to_tax, non_mono_to_register);
            }
            is_taxon_ptr->get_data().repr_in_solution_by = root_solution_node;
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
        n->get_data().repr_in_solution_by = root_solution_node;
        n->get_data().des_ids.clear();
    }
    for (auto n : postOrderInTaxNd) {
        n->get_data().des_ids.clear();
    }
    set<N *> visited;
    for (auto n : solnLeaves) {
        auto c = n;
        while (visited.count(c) == 0 && c != root_solution_node) {
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
    root_solution_node->get_data().des_ids.clear();
    root_solution_node->get_data().des_ids.insert(ott_id);
    // here we need to update the altered inc_sed taxon maps...
    // All the ones finished in this slice can be deleted from maps...
    for (auto is_taxon_ptr : inc_sed_that_are_proper_children) {
        unprune_stats.inc_sed_taxon_to_sampled_tips.at(is_taxon_ptr).clear();
        if (is_taxon_ptr->get_parent()) {
            is_taxon_ptr->detach_this_node();
        }
    }
    for (auto is_it = inc_sed_to_deal_with_by_level.rbegin(); is_it != inc_sed_to_deal_with_by_level.rend(); ++is_it) {
        for (auto is_taxon_ptr : is_it->second) {
            unprune_stats.inc_sed_taxon_to_sampled_tips.at(is_taxon_ptr).clear();
            if (is_taxon_ptr->get_parent()) {
                is_taxon_ptr->detach_this_node();
            }
        }
    }
    // all of the ones with MRCA deeper need to have the root_solution_node registered as the "sampled tip" (even though it really isn't a tip)
    for (auto is_taxon_ptr : inc_sed_mapping_deeper) {
        auto & samp_tip_set = unprune_stats.inc_sed_taxon_to_sampled_tips.at(is_taxon_ptr);
        auto & stsp = curr_slice_inc_sed_map[is_taxon_ptr];
        for (TaxSolnNdPair tsp : stsp) {
            // for the deeper parts of the tree the tips in this slice need to be annotated
            //    to reflect the fact that they have already been included in a taxonomic node
            Node_t * spike_nd = tsp.first;
            root_solution_node->get_data().des_ids.insert(spike_nd->get_ott_id());
            //flagIncSedAs
            while (true) {
                if (inc_sed_mapping_deeper.count(spike_nd) > 0) {
                    break;
                }
                spike_nd->get_data().repr_in_solution_by = rootTaxonNd;
                spike_nd = spike_nd->get_parent();
                assert(spike_nd != nullptr);
            }
            samp_tip_set.erase(tsp.second);
        }
        samp_tip_set.insert(root_solution_node);
    }
    unprune_stats.numInternalsAddedToBackbone += nodes_add_for_taxa.size();
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
        if (ott_to_tax.count(ist_id) == 0) {
            throw OTCError() << "Incertae sedis ID (" << ist_id << ") not in taxonomy.";
        }
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
}



// Taxa that are descendants of an incertae sedis taxon, will demand some special handling..
// so we find out which ones are not found in the solution before calling unprunePreppedInputs
// results are stored in unprune_stats. This relies on having been filled...
void findBrokenIncertaeSedisDescendants(Tree_t & taxonomy,
                                        Tree_t & solution,
                                        const map<long, childParPair> & ott_to_tax,
                                        UnpruneStats & unprune_stats) {
    const auto & to_tips_map = unprune_stats.inc_sed_taxon_to_sampled_tips;
    for (auto ttm_it: to_tips_map) {
        auto inc_sed_taxon = ttm_it.first;
        const auto & sampled = ttm_it.second;
        if (sampled.size() < 2) {
            continue;
        }
        const auto nonexc = inc_sed_taxon->get_data().nonexcluded_ids;
        const auto taxon_des = inc_sed_taxon->get_data().des_ids;
        assert(nonexc != nullptr);
        auto soln_mrca = find_mrca_via_traversing(sampled);
        for (auto t : iter_leaf_n(*soln_mrca)) {
            if (sampled.count(t) == 0 && 0 == nonexc->count(t->get_ott_id())) {
                unprune_stats.non_monophyletic_inc_sed[inc_sed_taxon] = soln_mrca;
                unprune_stats.non_monophyletic_inc_sed_to_desIds[inc_sed_taxon->get_ott_id()] = inc_sed_taxon->get_data().des_ids;
                break;
            }
        }
    }
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
