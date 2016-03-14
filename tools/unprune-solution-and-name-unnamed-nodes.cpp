#include <algorithm>
#include <unordered_set>
#include <iterator>
#include <sstream>
#include <stack>
#include <tuple>

#include "otc/otcli.h"
#include "otc/node_naming.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"

using namespace otc;
using std::vector;
using std::unique_ptr;
using std::string;
using std::map;

static std::string lostTaxaJSONFilename;

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
    private:
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
    assert(false);
}

// adds ottId to the desIds field of every node from firstNd down to ancAndLast (inclusive)
template <typename N>
void addToDesIdsForAnc(long ottId, N *firstNd, N *ancAndLast) {
    assert(firstNd);
    while (true) {
        firstNd->getData().desIds.insert(ottId);
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
void incorporateHigherTaxonNode(N* higherTaxonNd,
                                N* rootSolnNd,
                                std::set<N *> & nodesAddedForTaxa,
                                LostTaxonMap & ltm) {
    assert(higherTaxonNd);
    assert(rootSolnNd);
    const auto & taxDes = higherTaxonNd->getData().desIds;
    assert(taxDes.size() > 1);
    N * currSolnNd = rootSolnNd;
    N * tipmostMRCA = nullptr;
    // the root of the solution, must have all of the constituent taxa
    assert(isSubset(taxDes, rootSolnNd->getData().desIds));
    while (true) {
        const auto & solDes = currSolnNd->getData().desIds;
        if (solDes == taxDes) {
            if (higherTaxonNd->isOutDegreeOneNode()) {
                auto nt = introduceMonotypicParent(higherTaxonNd, currSolnNd);
                nodesAddedForTaxa.insert(nt);
            } else if (nodesAddedForTaxa.count(currSolnNd) > 0) {
                auto nt = addParentAndMoveUnsampledTaxChildren(higherTaxonNd, currSolnNd);
                nodesAddedForTaxa.insert(nt);
            } else {
                moveUnsampledChildren(higherTaxonNd, currSolnNd);
            }
            return;
        }
        assert(isProperSubset(taxDes, solDes));
        // look in this subtree for the higher taxon 
        tipmostMRCA = currSolnNd;
        N * nextSolnNd = nullptr;
        for (auto c : iter_child(*currSolnNd)) {
            const auto & cdesId = c->getData().desIds;
            if (!areDisjoint(cdesId, taxDes)) {
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
                    return;
                }
            }
        }
        assert(nextSolnNd != nullptr);
        currSolnNd = nextSolnNd;
    }
}

template <typename N>
void unpruneTaxaForSubtree(N *rootSolnNd,
                           map<long, N*> & ott2tax, 
                           map<long, N*> & ott2soln, 
                           LostTaxonMap & ltm) {
    assert(rootSolnNd);
    assert(rootSolnNd->hasOttId());
    assert(!rootSolnNd->isTip());
    const auto ottId = rootSolnNd->getOttId();
    assert(ott2soln.at(ottId) == rootSolnNd);
    N * rootTaxonNd = ott2tax.at(ottId);
    assert(!rootTaxonNd->isTip());
    // this will have the IDs for all of the include taxa for this slice of the tree
    auto & solnDesIds = rootSolnNd->getData().desIds;
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
    // If a tip of the solution is a higher taxon, then we should
    //  graft on the other tips here...
    for (auto l : solnLeaves) {
        if (l->isTip()) {
            assert(l->hasOttId());
            auto leafOttId = l->getOttId();
            auto taxonForLeaf = ott2tax.at(leafOttId);
            if (!taxonForLeaf->isTip()) {
                auto cvec = all_children(taxonForLeaf);
                for (auto c : cvec) {
                    taxonForLeaf->removeChild(c);
                    l->addChild(c);
                }
            }
        }
    }
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
        incorporateHigherTaxonNode(higherTaxonNd, rootSolnNd, nodesAddedForTaxa, ltm);
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
}


// Note that this function will steal some nodes from `taxonomy`
//  so, on output that will no longer be a valid tree. This should
//  be fine if both taxonomy and solution go out of scope together.
template<typename T>
LostTaxonMap unpruneTaxa(T & taxonomy, T & solution) {
    LostTaxonMap ltm;
    typedef typename T::node_type N;
    std::cerr << "Leaves:           solution = " << countLeaves(solution) << "   taxonomy = " << countLeaves(taxonomy) << std::endl;
    std::cerr << "Internal:         solution = " << n_internal(solution) << "   taxonomy = " << n_internal(taxonomy) << std::endl;
    std::cerr << "Internal splits:  solution = " << n_internal_out_degree_many(solution) << "   taxonomy = " << n_internal_out_degree_many(taxonomy) << std::endl;
    const auto out_degree_many1 = n_internal_out_degree_many(taxonomy);
    const auto n_internal_confirmed = n_internal_with_ott_id(solution);
    const auto n_internal_new = n_internal(solution) - n_internal_confirmed;
    // 1. First, remove nodes from the taxonomy that do not occur in the solution
    // 1a. Index solution nodes by OttId.
    map<long, N*> ott_to_tax;
    for (auto nd: iter_post(taxonomy)){
        if (nd->hasOttId()) {
            ott_to_tax[nd->getOttId()] = nd;
        } else {
            LOG(WARNING) << "  warning: node in taxonomy without an OTT ID.\n";
        }
    }
    map<long, N*> ott_to_sol;
    const auto snVec = all_internal_nodes_post(solution);
    // postorder walk over solution. Every time we find a taxon assigned to a taxon
    //  we augment the slice of the tree that is rooted at that node (and is the
    //  subtree that is cut at the deepest taxonomic node)
    for (auto nd: snVec){
        if (nd->isTip()) {
            if (!nd->hasOttId()) {
                throw OTCError() << "Tip "<< nd->getName() << " in solution lacks an OTT ID";
            }
            if (not ott_to_tax.count(nd->getOttId())) {
                throw OTCError() << "OttId "<< nd->getOttId() << " not in taxonomy!";
            }
        }
        auto p = nd->getParent();
        if (nd->hasOttId()){
            ott_to_sol[nd->getOttId()] = nd;
            if (!nd->isTip()) {
                // unpruneTaxaForSubtree will deal with this slice of the
                //tree
                unpruneTaxaForSubtree(nd, ott_to_tax, ott_to_sol, ltm);
                // we treat this named node as a tip for the next deeper slice
                nd->getData().desIds.clear();
            }
            nd->getData().desIds.insert(nd->getOttId());
            if (p) {
                auto & pd = p->getData();
                pd.desIds.insert(nd->getOttId());
            }
        } else if (p) {
            const auto & d = nd->getData().desIds;
            p->getData().desIds.insert(d.begin(), d.end());
        }
    }    
    return ltm;
}

bool handleLostTaxaJSON(OTCLI&, const std::string & arg) {
    lostTaxaJSONFilename = arg;
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
    
    auto & solution = *(trees.at(0));
    auto & taxonomy = *(trees.at(1));
    const auto lostTaxa = unpruneTaxa(taxonomy, solution);
    nameUnamedNodes(solution);
    writeTreeAsNewick(std::cout, solution);
    if (!lostTaxaJSONFilename.empty()) {
        writeLostTaxa(lostTaxaJSONStream, lostTaxa);
    }
    return 0;
}
