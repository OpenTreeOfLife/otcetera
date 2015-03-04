#include "otc/otcli.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef otc::RootedTree<typename Node_t::data_type, RTreeOttIDMapping<typename Node_t::data_type>> Tree_t;
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);
bool handleTabooOTTIdListFile(OTCLI & otCLI, const std::string &nextArg);

struct NodePairing {
    Node_t * scaffoldNode;
    Node_t * phyloNode;
    NodePairing(Node_t *taxo, Node_t *phylo)
        :scaffoldNode(taxo),
        phyloNode(phylo) {
    }
};
struct PathPairing {
    const Node_t * const phyloChild;
    const Node_t * const phyloParent;
    const Node_t * const scaffoldDes;
    const Node_t * const scaffoldAnc;

    PathPairing(const NodePairing & parent, const NodePairing & child)
        :phyloChild(child.phyloNode),
        phyloParent(parent.phyloNode),
        scaffoldDes(child.scaffoldNode),
        scaffoldAnc(parent.scaffoldNode) {
    }
};
struct AlignmentThreading {
    std::map<Tree_t *, std::set<NodePairing *> > nodeAlignments;
    std::map<Tree_t *, std::set<PathPairing *> > edgeBelowAlignments;
    std::map<Tree_t *, std::set<PathPairing *> > loopAlignments;
};

struct RemapToDeepestUnlistedState {
    std::unique_ptr<Tree_t> taxonomy;
    int numErrors;
    std::set<long> ottIds;
    std::set<const Node_t *> contestedNodes;
    std::set<long> tabooIds;
    std::list<std::unique_ptr<Tree_t> > inputTrees;
    std::list<NodePairing> nodePairings;
    std::list<PathPairing> pathPairings;
    std::map<const Node_t*, AlignmentThreading> taxoToAlignment;


    RemapToDeepestUnlistedState()
        :taxonomy(nullptr),
         numErrors(0) {
    }

    void summarize(const OTCLI &otCLI) {
        assert (taxonomy != nullptr);
        for (auto nd : contestedNodes) {
            otCLI.out << nd->getName() << '\n';
        }
    }

    bool processTaxonomyTree(OTCLI & otCLI) {
        ottIds = keys(taxonomy->getData().ottIdToNode);
        for (auto nd : iter_node(*taxonomy)) {
            taxoToAlignment.emplace(nd, AlignmentThreading{});
        }
        otCLI.getParsingRules().ottIdValidator = &ottIds;
        return true;
    }

    NodePairing * _addNodeMapping(Node_t *taxo, Node_t *nd, Tree_t *tree) {
        nodePairings.emplace_back(taxo, nd);
        auto ndPairPtr = &(*nodePairings.rbegin());
        auto & athreading = taxoToAlignment[taxo];
        athreading.nodeAlignments[tree].insert(ndPairPtr);
        return ndPairPtr;
    }
    PathPairing * _addPathMapping(NodePairing * parentPairing, NodePairing * childPairing, Tree_t *tree) {
        pathPairings.emplace_back(*parentPairing, *childPairing);
        auto pathPairPtr = &(*pathPairings.rbegin());
        // register a pointer to the path at each traversed...
        auto currTaxo = pathPairPtr->scaffoldDes;
        auto ancTaxo = pathPairPtr->scaffoldAnc;
        if (currTaxo != ancTaxo) {
            while (currTaxo != ancTaxo) {
                taxoToAlignment[currTaxo].edgeBelowAlignments[tree].insert(pathPairPtr);
                currTaxo = currTaxo->getParent();
                if (currTaxo == nullptr) {
                    break;
                }
            }
        } else {
            taxoToAlignment[currTaxo].loopAlignments[tree].insert(pathPairPtr);
        }
        return pathPairPtr;
    }
    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
        assert(taxonomy != nullptr);
        assert(tree != nullptr);
        std::map<Node_t *, NodePairing *> currTreeNodePairings;
        for (auto nd : iter_post(*tree)) {
            auto par = nd->getParent();
            if (par == nullptr) {
                continue;
            }
            NodePairing * ndPairPtr = nullptr;
            Node_t * taxoChild = nullptr;
            auto reuseNodePairingIt = currTreeNodePairings.find(nd);
            if (reuseNodePairingIt == currTreeNodePairings.end()) {
                auto ottId = nd->getOttId();
                taxoChild = taxonomy->getData().getNodeForOttId(ottId);
                ndPairPtr = _addNodeMapping(taxoChild, nd, tree.get());
                currTreeNodePairings[nd] = ndPairPtr;
            } else {
                ndPairPtr = reuseNodePairingIt->second;
                taxoChild = ndPairPtr->scaffoldNode;
            }
            NodePairing * parPairPtr = nullptr;
            auto prevAddedNodePairingIt = currTreeNodePairings.find(par);
            if (prevAddedNodePairingIt == currTreeNodePairings.end()) {
                const auto & parDesIds = par->getData().desIds;
                auto taxoPar = searchAncForMRCAOfDesIds(taxoChild, parDesIds);
                parPairPtr = _addNodeMapping(taxoPar, par, tree.get());
                currTreeNodePairings[par] = parPairPtr;
            } else {
                parPairPtr = prevAddedNodePairingIt->second;
            }
            _addPathMapping(parPairPtr, ndPairPtr, tree.get());
        }
        return true;
    }


    bool parseAndTabooOTTIdListFile(const std::string &fp) {
        auto t = parseListOfOttIds(fp);
        tabooIds.insert(begin(t), end(t));
        return true;
    }
};

inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
    RemapToDeepestUnlistedState * ctsp = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
    assert(ctsp != nullptr);
    assert(tree != nullptr);
    if (ctsp->taxonomy == nullptr) {
        ctsp->taxonomy = std::move(tree);
        return ctsp->processTaxonomyTree(otCLI);
    }
    return ctsp->processSourceTree(otCLI, std::move(tree));
}

bool handleTabooOTTIdListFile(OTCLI & otCLI, const std::string &nextArg) {
    RemapToDeepestUnlistedState * fusp = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
    assert(fusp != nullptr);
    assert(!nextArg.empty());
    return fusp->parseAndTabooOTTIdListFile(nextArg);
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otcdetectcontested",
                "takes at least 2 newick file paths: a full taxonomy tree, and some number of input trees. Writes the OTT IDs of clades in the taxonomy whose monophyly is questioned by at least one input",
                "taxonomy.tre inp1.tre inp2.tre");
    RemapToDeepestUnlistedState fus;
    otCLI.blob = static_cast<void *>(&fus);
    otCLI.addFlag('m',
                  "ARG=a file containing a list of taboo OTT ids.",
                  handleTabooOTTIdListFile,
                  true);
    auto rc = treeProcessingMain<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 2);
    if (rc == 0) {
        fus.summarize(otCLI);
        return fus.numErrors;
    }
    return rc;
}
