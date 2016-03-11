#include <algorithm>
#include <unordered_set>
#include <iterator>
#include <sstream>

#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_iter.h"
#include "otc/node_naming.h"
using namespace otc;
using std::vector;
using std::unique_ptr;
using std::string;

static std::string lostTaxaJSONFilename;

struct RTNodeSmallestChild {
    int smallestChild = 0;
};

using Tree_t = RootedTree<RTNodeSmallestChild, RTreeNoData>;

inline int smallestChild(const Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

inline int& smallestChild(Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

class LostTaxonLocation {

};
typedef std::map<int, LostTaxonLocation> LostTaxonMap;

template<typename T>
LostTaxonMap unpruneTaxa(T & taxonomy, T & solution);

void writeLostTaxa(std::ostream & out, const LostTaxonMap & ltm);

void writeLostTaxa(std::ostream & out, const LostTaxonMap & ltm) {
}
template<typename T>
LostTaxonMap unpruneTaxa(T & taxonomy, T & solution) {
    LostTaxonMap ltm;
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
