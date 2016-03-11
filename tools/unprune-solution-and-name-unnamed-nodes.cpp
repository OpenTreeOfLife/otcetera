#include <algorithm>
#include <unordered_set>
#include <iterator>
#include <sstream>

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

struct RTNodeSmallestChild {
    long smallestChild = 0;
};

using Tree_t = RootedTree<RTNodeSmallestChild, RTreeNoData>;

inline long smallestChild(const Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

inline long& smallestChild(Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

class LostTaxonLocation {
    public:
        void addAttachmentPoint(Tree_t::node_type *);
};
typedef std::map<int, LostTaxonLocation> LostTaxonMap;

template<typename T>
LostTaxonMap unpruneTaxa(T & taxonomy, T & solution);

void writeLostTaxa(std::ostream & out, const LostTaxonMap & ltm);

void writeLostTaxa(std::ostream & out, const LostTaxonMap & ltm) {
}

template <typename N>
void stealNodeAsChild(N * newPar, N *child) {

}


template <typename N>
N * mrcaForThisTaxon(N * parMRCASolnNd,
                     N * taxonNd,
                     map<long, N*> & ott2tax, 
                     map<long, N*> & ott2soln, 
                     LostTaxonMap & ltm) {

}

template <typename N>
void recursiveUnpruneTaxonForSubtree(N * nearestSolnAncWithTaxon,
                                     N * mrcaSolnNd,
                                     N * taxonNd,
                                     map<long, N*> & ott2tax, 
                                     map<long, N*> & ott2soln, 
                                     LostTaxonMap & ltm) {
    const auto ottId = taxonNd->getOttId();
    if (ott2soln.count(ottId) > 0) {
        return;
    }
    if (taxonNd->isTip()) {
        const auto taxPar = taxonNd->getParent();
        const auto taxParOttId = taxPar->getOttId();
        auto & ltl = ltm[taxParOttId];
        ltl.addAttachmentPoint(mrcaSolnNd);
        stealNodeAsChild(mrcaSolnNd, taxonNd);
    } else {
        N * nm = mrcaForThisTaxon(mrcaSolnNd, taxonNd, ott2tax, ott2soln, ltm);
        N * currTaxonChild = taxonNd->getFirstChild();
        assert(currTaxonChild != nullptr);
        N * nextTaxonChild = currTaxonChild->getNextSib();
        while (currTaxonChild) {
            recursiveUnpruneTaxonForSubtree(nearestSolnAncWithTaxon,
                                            nm,
                                            currTaxonChild,
                                            ott2tax,
                                            ott2soln, 
                                            ltm);
            currTaxonChild = nextTaxonChild;
            nextTaxonChild = currTaxonChild->getNextSib();
        }
    }
}

template <typename N>
void unpruneTaxaForSubtree(N *solnNd,
                           map<long, N*> & ott2tax, 
                           map<long, N*> & ott2soln, 
                           LostTaxonMap & ltm) {
    assert(solnNd);
    assert(solnNd->hasOttId());
    assert(!solnNd->isTip());
    const auto ottId = solnNd->getOttId();
    assert(ott2soln.at(ottId) == solnNd);
    N * taxonNd = ott2tax.at(ottId);
    assert(!taxonNd->isTip());
    N * currTaxonChild = taxonNd->getFirstChild();
    assert(currTaxonChild != nullptr);
    N * nextTaxonChild = currTaxonChild->getNextSib();
    while (currTaxonChild) {
        recursiveUnpruneTaxonForSubtree(solnNd,
                                        solnNd,
                                        currTaxonChild,
                                        ott2tax,
                                        ott2soln, 
                                        ltm);
        currTaxonChild = nextTaxonChild;
        nextTaxonChild = currTaxonChild->getNextSib();
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
    for (auto nd: iter_post(solution)){
        if (nd->hasOttId()){
            ott_to_sol[nd->getOttId()] = nd;
            if (!nd->isTip()) {
                unpruneTaxaForSubtree(nd, ott_to_tax, ott_to_sol, ltm);
            }
        }
        if (nd->isTip()) {
            if (!nd->hasOttId()) {
                throw OTCError() << "Tip "<< nd->getName() << " in solution lacks an OTT ID";
            }
            if (not ott_to_tax.count(nd->getOttId())) {
                throw OTCError() << "OttId "<< nd->getOttId() << " not in taxonomy!";
            }

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
