#include <tuple>
#include "otc/otcli.h"
#include "otc/debug.h"
#include "otc/greedy_forest.h"
#include "otc/embedding.h"
using namespace otc;
std::unique_ptr<TreeMappedWithSplits> cloneTree(const TreeMappedWithSplits &);

template<typename T, typename U>
void updateAncestralPathOttIdSet(T * nd, const OttIdSet & oldEls, const OttIdSet newEls, std::map<const T *, NodeThreading<T, U> > & m);
template<typename T, typename U>
inline void updateAncestralPathOttIdSet(T * nd, const OttIdSet & oldEls, const OttIdSet newEls, std::map<const T *, NodeThreading<T, U> > & m) {
    for (auto anc : iter_anc(*nd)) {
         auto & ant = m.at(anc);
         ant.updateAllPathsOttIdSets(oldEls, newEls);
    }
}

//currently not copying names
std::unique_ptr<TreeMappedWithSplits> cloneTree(const TreeMappedWithSplits &tree) {
    TreeMappedWithSplits * rawTreePtr = new TreeMappedWithSplits();
    try {
        NodeWithSplits * newRoot = rawTreePtr->createRoot();
        auto r = tree.getRoot();
        assert(r->hasOttId());
        newRoot->setOttId(r->getOttId());
        std::map<const NodeWithSplits *, NodeWithSplits *> templateToNew;
        templateToNew[r]= newRoot;
        std::map<long, NodeWithSplits *> & newMap = rawTreePtr->getData().ottIdToNode;
        rawTreePtr->getData().desIdSetsContainInternals = tree.getData().desIdSetsContainInternals;
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
            }
            nn->getData().desIds = nd->getData().desIds;
        }
    } catch (...) {
        delete rawTreePtr;
        throw;
    }
    return std::unique_ptr<TreeMappedWithSplits>(rawTreePtr);
}



using NodePairingWithSplits = NodePairing<NodeWithSplits, NodeWithSplits>;
using PathPairingWithSplits = PathPairing<NodeWithSplits, NodeWithSplits>;
using NodeThreadingWithSplits = NodeThreading<NodeWithSplits, NodeWithSplits>;

class ThreadedTree {
    protected:
    std::list<NodePairingWithSplits> nodePairings;
    std::list<PathPairingWithSplits> pathPairings;
    std::map<const NodeWithSplits *, NodeThreadingWithSplits> taxoToAlignment;
    public:

    NodePairingWithSplits * _addNodeMapping(NodeWithSplits *taxo, NodeWithSplits *nd, std::size_t treeIndex) {
        assert(taxo != nullptr);
        assert(nd != nullptr);
        nodePairings.emplace_back(NodePairingWithSplits(taxo, nd));
        auto ndPairPtr = &(*nodePairings.rbegin());
        auto & athreading = taxoToAlignment[taxo];
        athreading.nodeAlignments[treeIndex].insert(ndPairPtr);
        return ndPairPtr;
    }
    PathPairingWithSplits * _addPathMapping(NodePairingWithSplits * parentPairing,
                                            NodePairingWithSplits * childPairing,
                                            std::size_t treeIndex) {
        pathPairings.emplace_back(*parentPairing, *childPairing);
        auto pathPairPtr = &(*pathPairings.rbegin());
        // register a pointer to the path at each traversed...
        auto currTaxo = pathPairPtr->scaffoldDes;
        auto ancTaxo = pathPairPtr->scaffoldAnc;
        if (currTaxo != ancTaxo) {
            while (currTaxo != ancTaxo) {
                taxoToAlignment[currTaxo].edgeBelowAlignments[treeIndex].insert(pathPairPtr);
                currTaxo = currTaxo->getParent();
                if (currTaxo == nullptr) {
                    break;
                }
            }
        } else {
            taxoToAlignment[currTaxo].loopAlignments[treeIndex].insert(pathPairPtr);
        }
        return pathPairPtr;
    }
    void threadNewTree(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex) {
        // do threading
        std::map<NodeWithSplits *, NodePairingWithSplits *> currTreeNodePairings;
        std::set<NodePairingWithSplits *> tipPairings;
        for (auto nd : iter_post(tree)) {
            auto par = nd->getParent();
            if (par == nullptr) {
                continue;
            }
            NodePairingWithSplits * ndPairPtr = nullptr;
            NodeWithSplits * taxoDes = nullptr;
            if (nd->isTip()) {
                assert(currTreeNodePairings.find(nd) == currTreeNodePairings.end()); // TMP, Remove this to save time?
                assert(nd->hasOttId());
                auto ottId = nd->getOttId();
                taxoDes = scaffoldTree.getData().getNodeForOttId(ottId);
                assert(taxoDes != nullptr);
                ndPairPtr = _addNodeMapping(taxoDes, nd, treeIndex);
                for (auto former : tipPairings) {
                    if (areLinearlyRelated(taxoDes, former->scaffoldNode)) {
                        std::string m = "Repeated or nested OTT ID in tip mapping of an input tree: \"";
                        m += nd->getName();
                        m += "\" and \"";
                        m += former->phyloNode->getName();
                        m += "\" found.";
                        throw OTCError(m);
                    }
                }
                tipPairings.insert(ndPairPtr);
                currTreeNodePairings[nd] = ndPairPtr;
            } else {
                auto reuseNodePairingIt = currTreeNodePairings.find(nd);
                assert(reuseNodePairingIt != currTreeNodePairings.end());
                ndPairPtr = reuseNodePairingIt->second;
                taxoDes = ndPairPtr->scaffoldNode;
                assert(taxoDes != nullptr);
            }
            NodePairingWithSplits * parPairPtr = nullptr;
            auto prevAddedNodePairingIt = currTreeNodePairings.find(par);
            if (prevAddedNodePairingIt == currTreeNodePairings.end()) {
                const auto & parDesIds = par->getData().desIds;
                auto taxoAnc = searchAncForMRCAOfDesIds(taxoDes, parDesIds);
                assert(taxoAnc != nullptr);
                parPairPtr = _addNodeMapping(taxoAnc, par, treeIndex);
                currTreeNodePairings[par] = parPairPtr;
            } else {
                parPairPtr = prevAddedNodePairingIt->second;
            }
            _addPathMapping(parPairPtr, ndPairPtr, treeIndex);
        }
    }

    void threadTaxonomyClone(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex) {
        // do threading
        std::map<NodeWithSplits *, NodePairingWithSplits *> currTreeNodePairings;
        for (auto nd : iter_post(tree)) {
            auto par = nd->getParent();
            if (par == nullptr) {
                continue;
            }
            NodePairingWithSplits * ndPairPtr = nullptr;
            NodeWithSplits * taxoDes = nullptr;
            if (nd->isTip()) {
                assert(currTreeNodePairings.find(nd) == currTreeNodePairings.end()); // TMP, Remove this to save time?
                assert(nd->hasOttId());
                auto ottId = nd->getOttId();
                taxoDes = scaffoldTree.getData().getNodeForOttId(ottId);
                assert(taxoDes != nullptr);
                ndPairPtr = _addNodeMapping(taxoDes, nd, treeIndex);
                currTreeNodePairings[nd] = ndPairPtr;
            } else {
                auto reuseNodePairingIt = currTreeNodePairings.find(nd);
                assert(reuseNodePairingIt != currTreeNodePairings.end());
                ndPairPtr = reuseNodePairingIt->second;
                taxoDes = ndPairPtr->scaffoldNode;
                assert(taxoDes != nullptr);
            }
            NodePairingWithSplits * parPairPtr = nullptr;
            auto prevAddedNodePairingIt = currTreeNodePairings.find(par);
            if (prevAddedNodePairingIt == currTreeNodePairings.end()) {
                auto pottId = par->getOttId(); // since it is a taxonomy, it will have internal node labels
                auto taxoAnc = scaffoldTree.getData().getNodeForOttId(pottId);
                assert(taxoAnc != nullptr);
                parPairPtr = _addNodeMapping(taxoAnc, par, treeIndex);
                currTreeNodePairings[par] = parPairPtr;
            } else {
                parPairPtr = prevAddedNodePairingIt->second;
            }
            _addPathMapping(parPairPtr, ndPairPtr, treeIndex);
        }
    }

    void writeDOTExport(std::ostream & out,
                           const NodeThreading<NodeWithSplits, NodeWithSplits> & thr,
                           const NodeWithSplits * nd,
                           const std::vector<TreeMappedWithSplits *> & ) const {
        writeNewick(out, nd);
        out << "nd.ottID = " << nd->getOttId() << " --> " << (nd->getParent() ? nd->getParent()->getOttId() : 0L) << "\n";
        out << "  getTotalNumNodeMappings = " << thr.getTotalNumNodeMappings() << "\n";
        out << "  getTotalNumLoops = " << thr.getTotalNumLoops() << "\n";
        out << "  getTotalNumEdgeBelowTraversals = " << thr.getTotalNumEdgeBelowTraversals() << "\n";
        out << "  isContested = " << thr.isContested() << "\n";
    }

};


class RemapToDeepestUnlistedState
    : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits>,
    public ThreadedTree {
    public:
    int numErrors;
    OttIdSet tabooIds;
    std::map<std::unique_ptr<TreeMappedWithSplits>, std::size_t> inputTreesToIndex;
    std::vector<TreeMappedWithSplits *> treePtrByIndex;
    bool doReportAllContested;
    bool doConstructSupertree;
    std::list<long> idsListToReportOn;
    std::list<long> idListForDotExport;
    TreeMappedWithSplits * taxonomyAsSource;

    void resolveOrCollapse(NodeWithSplits * scaffoldNd, SupertreeContextWithSplits & sc) {
        auto & thr = taxoToAlignment[scaffoldNd];
        if (thr.isContested()) {
            if (thr.highRankingTreesPreserveMonophyly(sc.numTrees)) {
                thr.resolveGivenContestedMonophyly(*scaffoldNd, sc);
            } else {
                thr.constructPhyloGraphAndCollapseIfNecessary(*scaffoldNd, sc);
            }
        } else {
            thr.resolveGivenUncontestedMonophyly(*scaffoldNd, sc);
        }
    }
    void constructSupertree() {
        const auto numTrees = treePtrByIndex.size();
        TreeMappedWithSplits * tax = taxonomy.get();
        SupertreeContextWithSplits sc{numTrees, taxoToAlignment, *tax};
        LOG(DEBUG) << "Before supertree "; writeTreeAsNewick(std::cerr, *taxonomy); std::cerr << '\n';
        for (auto nd : iter_post_internal(*taxonomy)) {
            if (nd == taxonomy->getRoot()) resolveOrCollapse(nd, sc);
            LOG(DEBUG) << "After handling " << nd->getOttId(); writeTreeAsNewick(std::cerr, *taxonomy); std::cerr << '\n';
        }
        
    }


    virtual ~RemapToDeepestUnlistedState(){}
    RemapToDeepestUnlistedState()
        :TaxonomyDependentTreeProcessor<TreeMappedWithSplits>(),
         numErrors(0),
         doReportAllContested(false),
         doConstructSupertree(false),
         taxonomyAsSource(nullptr) {
    }

    void reportAllConflicting(std::ostream & out, bool verbose) {
        std::map<std::size_t, unsigned long> nodeMappingDegree;
        std::map<std::size_t, unsigned long> passThroughDegree;
        std::map<std::size_t, unsigned long> loopDegree;
        unsigned long totalContested = 0;
        unsigned long redundContested = 0;
        unsigned long totalNumNodes = 0;
        for (auto nd : iter_node_internal(*taxonomy)) {
            const auto & thr = taxoToAlignment[nd];
            nodeMappingDegree[thr.getTotalNumNodeMappings()] += 1;
            passThroughDegree[thr.getTotalNumEdgeBelowTraversals()] += 1;
            loopDegree[thr.getTotalNumLoops()] += 1;
            totalNumNodes += 1;
            std::vector<NodeWithSplits *> aliasedBy = getNodesAliasedBy(nd, *taxonomy);
            if (thr.reportIfContested(out, nd, treePtrByIndex, aliasedBy, verbose)) {
                totalContested += 1;
                if (nd->getOutDegree() == 1) {
                    redundContested += 1;
                }
            }
        }
        unsigned long m = std::max(loopDegree.rbegin()->first, passThroughDegree.rbegin()->first);
        m = std::max(m, nodeMappingDegree.rbegin()->first);
        out << "Degree\tNodeMaps\tEdgeMaps\tLoops\n";
        for (unsigned long i = 0 ; i <= m; ++i) {
            out << i << '\t' << nodeMappingDegree[i]<< '\t' << passThroughDegree[i] << '\t' << loopDegree[i]<< '\n';
        }
        out << totalNumNodes << " internals\n" << totalContested << " contested\n" << (totalNumNodes - totalContested) << " uncontested\n";
        out << redundContested << " monotypic contested\n";
    }
    
    bool summarize(const OTCLI &otCLI) override {
        if (doConstructSupertree) {
            cloneTaxonomyAsASourceTree();
            constructSupertree();
            writeTreeAsNewick(otCLI.out, *taxonomy);
            otCLI.out << '\n';
        }
        std::ostream & out{otCLI.out};
        assert (taxonomy != nullptr);
        if (doReportAllContested) {
            reportAllConflicting(out, otCLI.verbose);
        } else {
            for (auto tr : idsListToReportOn) {
                auto nd = taxonomy->getData().getNodeForOttId(tr);
                if (nd == nullptr) {
                    throw OTCError(std::string("Unrecognized OTT ID in list of OTT IDs to report on: ") + std::to_string(tr));
                }
                const auto & thr = taxoToAlignment[nd];
                std::vector<NodeWithSplits *> aliasedBy = getNodesAliasedBy(nd, *taxonomy);
                thr.reportIfContested(out, nd, treePtrByIndex, aliasedBy, otCLI.verbose);
            }
        }
        for (auto tr : idListForDotExport) {
            auto nd = taxonomy->getData().getNodeForOttId(tr);
            if (nd == nullptr) {
                throw OTCError(std::string("Unrecognized OTT ID in list of OTT IDs to export to DOT: ") + std::to_string(tr));
            }
            //const auto & thr = taxoToAlignment[nd];
            //writeDOTExport(out, thr, nd, treePtrByIndex);
            for (auto n : iter_pre_n_const(nd)) {
                const auto & thr = taxoToAlignment[n];
                writeDOTExport(out, thr, n, treePtrByIndex);
            }
        }
        return true;
    }

    bool processTaxonomyTree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::processTaxonomyTree(otCLI);
        checkTreeInvariants(*taxonomy);
        suppressMonotypicTaxaPreserveDeepestDangle(*taxonomy);
        checkTreeInvariants(*taxonomy);
        for (auto nd : iter_node(*taxonomy)) {
            taxoToAlignment.emplace(nd, NodeThreadingWithSplits{});
        }
        return true;
    }

    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> treeup) override {
        assert(treeup != nullptr);
        assert(taxonomy != nullptr);
        // Store the tree pointer with a map to its index, and an alias for fast index->tree.
        std::size_t treeIndex = inputTreesToIndex.size();
        assert(treeIndex == treePtrByIndex.size());
        TreeMappedWithSplits * raw = treeup.get();
        inputTreesToIndex[std::move(treeup)] = treeIndex;
        treePtrByIndex.push_back(raw);
        // Store the tree's filename
        raw->setName(otCLI.currentFilename);
        threadNewTree(*taxonomy, *raw, treeIndex);
        otCLI.out << "# pathPairings = " << pathPairings.size() << '\n';
        return true;
    }

    bool cloneTaxonomyAsASourceTree() {
        assert(taxonomy != nullptr);
        assert(taxonomyAsSource == nullptr);
        std::unique_ptr<TreeMappedWithSplits> tree = std::move(cloneTree(*taxonomy));
        taxonomyAsSource = tree.get();
        std::size_t treeIndex = inputTreesToIndex.size();
        TreeMappedWithSplits * raw = tree.get();
        inputTreesToIndex[std::move(tree)] = treeIndex;
        treePtrByIndex.push_back(taxonomyAsSource);
        // Store the tree's filename
        raw->setName("TAXONOMY");
        threadTaxonomyClone(*taxonomy, *taxonomyAsSource, treeIndex);
        return true;
    }
};

bool handleReportAllFlag(OTCLI & otCLI, const std::string &);
bool handleReportOnNodesFlag(OTCLI & otCLI, const std::string &);
bool handleDotNodesFlag(OTCLI & otCLI, const std::string &narg);
bool handleSuperTreeFlag(OTCLI & otCLI, const std::string &narg);

bool handleReportAllFlag(OTCLI & otCLI, const std::string &) {
    RemapToDeepestUnlistedState * proc = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->doReportAllContested = true;
    return true;
}
bool handleSuperTreeFlag(OTCLI & otCLI, const std::string &) {
    RemapToDeepestUnlistedState * proc = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->doConstructSupertree = true;
    return true;
}

bool handleReportOnNodesFlag(OTCLI & otCLI, const std::string &narg) {
    RemapToDeepestUnlistedState * proc = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a list of IDs after the -b argument.");
    }
    auto rs = split_string(narg, ',');
    for (auto word : rs) {
        auto ottId = ottIDFromName(word);
        if (ottId < 0) {
            throw OTCError(std::string("Expecting a list of IDs after the -b argument. Offending word: ") + word);
        }
        proc->idsListToReportOn.push_back(ottId);
    }
    return true;
}

bool handleDotNodesFlag(OTCLI & otCLI, const std::string &narg) {
    RemapToDeepestUnlistedState * proc = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a list of IDs after the -d argument.");
    }
    auto rs = split_string(narg, ',');
    for (auto word : rs) {
        auto ottId = ottIDFromName(word);
        if (ottId < 0) {
            throw OTCError(std::string("Expecting a list of IDs after the -d argument. Offending word: ") + word);
        }
        proc->idListForDotExport.push_back(ottId);
    }
    return true;
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otcscaffoldedsupertree",
                "takes at least 2 newick file paths: a full taxonomy tree, and some number of input trees. Crashes or emits bogus output.",
                "taxonomy.tre inp1.tre inp2.tre");
    RemapToDeepestUnlistedState proc;
    otCLI.addFlag('a',
                  "Write a report of all contested nodes",
                  handleReportAllFlag,
                  false);
    otCLI.addFlag('s',
                  "Compute a supertree",
                  handleSuperTreeFlag,
                  false);
    otCLI.addFlag('b',
                  "IDLIST should be a list of OTT IDs. A status report will be generated for those nodes",
                  handleReportOnNodesFlag,
                  true);
    otCLI.addFlag('d',
                  "IDLIST should be a list of OTT IDs. A DOT file of the nodes will be generated ",
                  handleDotNodesFlag,
                  true);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}

