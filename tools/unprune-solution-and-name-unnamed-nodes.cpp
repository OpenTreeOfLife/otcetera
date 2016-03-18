#include <algorithm>
#include <unordered_set>
#include <iterator>
#include <sstream>
#include <stack>
#include <tuple>
#include "json.hpp"
#include "otc/otcli.h"
#include "otc/node_naming.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"

using namespace otc;
using json = nlohmann::json;
using std::vector;
using std::unique_ptr;
using std::string;
using std::map;

static std::string lostTaxaJSONFilename;
static std::string statsJSONFilename;

struct RTNodePartialDesSet {
    std::set<long> desIds;
    long smallestChild  = 0;
};

using Tree_t = RootedTree<RTNodePartialDesSet, RTreeNoData>;

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
LostTaxonMap unpruneTaxa(T & taxonomy, T & solution);

void writeLostTaxa(std::ostream & out, const LostTaxonMap & ltm);

void writeLostTaxa(std::ostream & out, const LostTaxonMap & ltm) {
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
    out << document.dump(4) << std::endl;
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
        firstNd = firstNd->getParent();
    }
}

template <typename N>
std::vector<N *> higherTaxPreOrderBelowBoundaries(N * root,
                                                  const std::set<N *> & boundaries) {
    std::set<N *> seen;
    std::vector<N *> r;
    N * curr = root;
    assert(!curr->getData().desIds.empty());
    std::stack<N *> toDealWith;
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

template <typename N>
void registerAttachmentPoints(const N * taxon,
                              N * mrca,
                              LostTaxonLocation & ltl) {
    assert(taxon);
    assert(mrca);
    std::set<N *> visited;
    std::stack<N *> toDealWith;
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
N * bisect_branch_with_new_child_copy_info(N * curr, N * infoSource) {
    N * nt = bisect_branch_with_new_child(curr);
    nt->setName(infoSource->getName());
    nt->setOttId(infoSource->getOttId());
    nt->getData().desIds = curr->getData().desIds;
    return nt;
}

// returns the number of nodes on the backbone that are merged with this higher taxon
//      because the taxon was compatible (the same as) the node - will be either 0 or 1.
template <typename N>
std::size_t incorporateHigherTaxonNode(N* higherTaxonNd,
                                N* rootSolnNd,
                                std::set<N *> & nodesAddedForTaxa,
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
                    nt = bisect_branch_with_new_child_copy_info(currSolnNd, higherTaxonNd);
                } else {
                    LOG(DEBUG) << "introduceMonotypicParent";
                    nt = introduceMonotypicParent(higherTaxonNd, currSolnNd);
                }
            } else if (currSolnNd->hasOttId() 
                      && currSolnNd->getOttId() != higherTaxonNd->getOttId()) {
                if (currSolnNd == rootSolnNd) {
                    LOG(DEBUG) << "Adding forking child to make parent monotypic";
                    nt = bisect_branch_with_new_child_copy_info(currSolnNd, higherTaxonNd);
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

using NumAddedTuple = std::tuple<std::size_t, std::size_t, std::size_t>;
// Moves unsampled taxa from the taxanomy (by finding them in `ott2tax`) to a slice of the
//  solution tree that is rooted at `rootSolnNd`.
// Uses the `rootSolnNd->getData().desIds` to denote the set of leaves of this slice of the
//  tree.
// If the taxon conflicts with the solution, then the unsampled taxa are attached at the
//  MRCA of the taxon in the solution, and an entry is added mapping the taxon's OTT Id
//  to a LostTaxonLocation object that explains where the elements of the taxon are located.
// returns a pair of numbers:
//      0 is number of new nodes added along the soln backbone (this does not include)
//          the internal nodes transferred from the taxonomy to the solution.
//      1 is the number of internal nodes that were merged with an existing node on the
//          backbone
//      2 is the number of solution leaves that were expanded to internals
template <typename N>
NumAddedTuple unpruneTaxaForSubtree(N *rootSolnNd,
                           map<long, N*> & ott2tax, 
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
    std::set<N *> solnLeaves;
    std::set<N *> taxaLeaves;
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
    std::size_t numExpanded = 0;
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
    std::size_t numMerged = numExpanded;
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
    std::set<N *> nodesAddedForTaxa;
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
    std::set<N *> visited;
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


// Note that this function will steal some nodes from `taxonomy`
//  so, on output that will no longer be a valid tree. This should
//  be fine if both taxonomy and solution go out of scope together.
template<typename T>
LostTaxonMap unpruneTaxa(T & taxonomy, T & solution, std::ostream * statsStreamPtr) {
    LostTaxonMap ltm;
    typedef typename T::node_type N;
    const auto out_degree_many1 = n_internal_out_degree_many(taxonomy);
    std::size_t startNumNamedInternalsInSoln = 0;
    std::size_t numUnnamedInternalsInSoln = 0;
    std::size_t numTaxaLeaves = 0;
    std::size_t numTaxaForkingInternals = 0;
    // 1. First, remove nodes from the taxonomy that do not occur in the solution
    // 1a. Index solution nodes by OttId.
    map<long, N*> ott_to_tax;
    OttIdSet monotypicOttIds;
    for (auto nd: iter_post(taxonomy)){
        if (nd->hasOttId()) {
            ott_to_tax[nd->getOttId()] = nd;
            if (nd->isOutDegreeOneNode()) {
                monotypicOttIds.insert(nd->getOttId());
            } else if (nd->isTip()) {
                numTaxaLeaves += 1;
            } else {
                numTaxaForkingInternals += 1;
            }
        } else {
            LOG(WARNING) << "  warning: node in taxonomy without an OTT ID.\n";
        }
    }
    const auto numTaxaMonotypicInternals = monotypicOttIds.size();
    map<long, N*> ott_to_sol;
    const auto snVec = all_nodes(solution);
    // postorder walk over solution. Every time we find a taxon assigned to a taxon
    //  we augment the slice of the tree that is rooted at that node (and is the
    //  subtree that is cut at the deepest taxonomic node)
    std::size_t startNumSolnLeaves = 0;
    std::size_t startNumSolnMonotypicInternals = 0;
    std::size_t startNumSolnForkingInternals = 0;
    std::size_t numInternalsAddedToBackbone = 0;
    std::size_t numInternalsMergedOnBackbone = 0;
    std::size_t numLeavesExpanded = 0;
    for (auto nd: snVec){
        if (nd->isTip()) {
            startNumSolnLeaves += 1;
            if (!nd->hasOttId()) {
                throw OTCError() << "Tip "<< nd->getName() << " in solution lacks an OTT ID";
            }
            if (not ott_to_tax.count(nd->getOttId())) {
                throw OTCError() << "OttId "<< nd->getOttId() << " not in taxonomy!";
            }
        } else {
            if (nd->isOutDegreeOneNode()) {
                startNumSolnMonotypicInternals += 1;
            } else {
                startNumSolnForkingInternals += 1;
            }
        }
        auto p = nd->getParent();
        if (nd->hasOttId()){
            ott_to_sol[nd->getOttId()] = nd;
            if (!nd->isTip()) {
                startNumNamedInternalsInSoln += 1;
                // unpruneTaxaForSubtree will deal with this slice of the
                //tree
                const auto x = unpruneTaxaForSubtree(nd, ott_to_tax, ott_to_sol, ltm);
                numInternalsAddedToBackbone += std::get<0>(x);
                numInternalsMergedOnBackbone += std::get<1>(x);
                numLeavesExpanded += std::get<2>(x);
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
                numUnnamedInternalsInSoln += 1;
            }
            if (p) {
                const auto & d = nd->getData().desIds;
                p->getData().desIds.insert(d.begin(), d.end());
            }
        }
    }
    // This is similar to, but different from, the number of non-monotypic nodes reject.
    // That is because the rejected nodes are marked as monotypic if they have no ANCESTRAL children.
    const auto numTaxaRejected = ltm.size();
    std::size_t numMonotypicTaxaRejected = 0U;
    for (auto monoOtt : monotypicOttIds) {
        if (ltm.count(monoOtt) > 0) {
            numMonotypicTaxaRejected += 1;
        }
    }
    const auto numForkingTaxaRejected = numTaxaRejected - numMonotypicTaxaRejected;
    const auto startNumSolnInternals = startNumSolnMonotypicInternals + startNumSolnForkingInternals;
    const auto numTaxaInternals = numTaxaMonotypicInternals + numTaxaForkingInternals;
    const auto numInternalsAddedByUnpruning = numTaxaInternals - numTaxaRejected - numInternalsMergedOnBackbone;
    const auto numLeavesAdded = numTaxaLeaves - startNumSolnLeaves;
    std::cerr << "Leaves:           solution = " << startNumSolnLeaves << "   taxonomy = " << numTaxaLeaves << std::endl;
    std::cerr << "Internal:         solution = " << startNumSolnInternals << "   taxonomy = " << numTaxaInternals << std::endl;
    std::cerr << "Internal splits:  solution = " << startNumSolnForkingInternals << "   taxonomy = " << numTaxaForkingInternals << std::endl;
    std::cerr << "Taxa rejected: = " << numTaxaRejected << std::endl;
    std::cerr << "Taxonomy splits: #rejected by phylo inputs = " << numForkingTaxaRejected << std::endl;
    std::cerr << "Solution splits: input # with OTT Ids       = " << startNumNamedInternalsInSoln << std::endl;
    if (statsStreamPtr) {
        json document;
        json input;
        json output;
        input["num_solution_leaves"] = startNumSolnLeaves;
        input["num_taxonomy_leaves"] = numTaxaLeaves;
        input["num_solution_internals"] = startNumSolnInternals;
        input["num_taxonomy_internals"] = numTaxaInternals;
        input["num_solution_splits"] = startNumSolnForkingInternals;
        input["num_taxonomy_splits"] = numTaxaForkingInternals;
        output["num_taxa_rejected"] = numTaxaRejected;
        output["num_taxonomy_splits_rejected"] = numForkingTaxaRejected;
        output["num_taxonomy_internals_merged"] = numInternalsMergedOnBackbone;
        output["num_solution_leaves_expanded"] = numLeavesExpanded;
        output["num_leaves_added"] = numLeavesAdded;
        document["input"] = input;
        document["output"] = output;
        *statsStreamPtr << document.dump(4) << std::endl;
    }
    return ltm;
}

bool handleLostTaxaJSON(OTCLI&, const std::string & arg) {
    lostTaxaJSONFilename = arg;
    return true;
}

bool handleStatsJSON(OTCLI&, const std::string & arg) {
    statsJSONFilename = arg;
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
   
    auto & solution = *(trees.at(0));
    auto & taxonomy = *(trees.at(1));
    const auto lostTaxa = unpruneTaxa(taxonomy, solution, statsStreamPtr);
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
