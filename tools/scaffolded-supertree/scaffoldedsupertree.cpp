#include <tuple>
#include "otc/otcli.h"
#include "otc/debug.h"
#include "otc/greedy_forest.h"
#include "otc/embedding.h"
#include "otc/embedded_tree.h"
using namespace otc;
std::unique_ptr<TreeMappedWithSplits> cloneTree(const TreeMappedWithSplits &);

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

enum SuperTreeDOTStep {
    BEFORE_ND_WO_TAXO,
    BEFORE_ND_W_TAXO,
    AFTER_ND_WO_TAXO,
    AFTER_ND_W_TAXO
};

class ScaffoldedSupertree
    : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits>,
    public EmbeddedTree {
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
    int currDotFileIndex;
    bool debuggingOutput;
    bool emitScaffoldDotFiles;

    void writeEmbeddingDOT(SuperTreeDOTStep sts, const NodeWithSplits * nd, const NodeWithSplits * actionNd) {
        std::string fn = "ScaffSuperTree_num";
        fn += std::to_string(currDotFileIndex++);
        fn += "_";
        const bool includeLastTree = (sts == BEFORE_ND_W_TAXO || sts == AFTER_ND_W_TAXO);
        if (sts == BEFORE_ND_WO_TAXO || sts == BEFORE_ND_W_TAXO) {
            fn += "Nd";
            fn += std::to_string(nd->getOttId());
            fn += "Before";
        } else if (sts == AFTER_ND_WO_TAXO || sts == AFTER_ND_W_TAXO) {
            fn += "Nd";
            fn += std::to_string(nd->getOttId());
            fn += "After";
        }
        if (actionNd != nullptr) {
            fn += "Nd";
            fn += std::to_string(actionNd->getOttId());
        }
        if (includeLastTree) {
            fn += "WTax";
        }
        fn += ".dot";
        std::ofstream out;
        out.open(fn);
        const auto & thr = _getEmdeddingForNode(nd);
        writeDOTExport(out, thr, nd, treePtrByIndex, true, includeLastTree);
    }
    void writeNumberedDOT(NodeWithSplits * nd, bool entireSubtree, bool includeLastTree) {
        std::string fn = "ScaffSuperTree" + std::to_string(currDotFileIndex++) + ".dot";
        LOG(DEBUG) << "writing DOT file \"" << fn << "\"";
        std::ofstream out;
        out.open(fn);
        const auto & thr = _getEmdeddingForNode(nd);
        writeDOTExport(out, thr, nd, treePtrByIndex, entireSubtree, includeLastTree);
        LOG(DEBUG) << "finished DOT file \"" << fn << "\"";
    }

    void resolveOrCollapse(NodeWithSplits * scaffoldNd, SupertreeContextWithSplits & sc) {
        assert(!scaffoldNd->isTip());
        auto & thr = _getEmdeddingForNode(scaffoldNd);
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
        SupertreeContextWithSplits sc{numTrees, taxoToEmbedding, *tax};
        if (debuggingOutput) {
            if (emitScaffoldDotFiles) {
                writeEmbeddingDOT(BEFORE_ND_W_TAXO, taxonomy->getRoot(), nullptr);
                writeEmbeddingDOT(BEFORE_ND_WO_TAXO, taxonomy->getRoot(), nullptr);
            } else {
                LOG(DEBUG) << "Beginning construction of the supertree from the embedded tree.";
            }
        }
        std::list<NodeWithSplits * > postOrder;
        for (auto nd : iter_post_internal(*taxonomy)) {
            assert(!nd->isTip());
            postOrder.push_back(nd);
        }
        for (auto nd : postOrder) {
            auto p = nd->getParent();
            if (debuggingOutput) {
                if (emitScaffoldDotFiles) {
                    writeEmbeddingDOT(BEFORE_ND_W_TAXO, nd, nd);
                    writeEmbeddingDOT(BEFORE_ND_WO_TAXO, nd, nd);
                } else {
                     
                }
            }
            LOG(INFO) << "Calling resolveOrCollapse for OTT" << nd->getOttId();
            resolveOrCollapse(nd, sc);
            if (debuggingOutput) {
                if (emitScaffoldDotFiles) {
                    writeEmbeddingDOT(AFTER_ND_W_TAXO, nd, nd);
                    writeEmbeddingDOT(AFTER_ND_WO_TAXO, nd, nd);
                    writeEmbeddingDOT(AFTER_ND_W_TAXO, taxonomy->getRoot(), nd);
                    writeEmbeddingDOT(AFTER_ND_WO_TAXO, taxonomy->getRoot(), nd);
                } else {
                    LOG(DEBUG) << "Completed resolveOrCollapse call for OTT" << nd->getOttId();
                }
            }
            if (p != nullptr) {
                checkAllNodePointersIter(p);
            }
            bool before = true;
            for (auto u : postOrder) {
                if (u == nd) {
                    before = false;
                } else if ((!before) && u->isTip()) {
                    LOG(ERROR) << "Node for OTT " << u->getOttId() << " has become a tip after processing OTT" << nd->getOttId();
                    assert(false);
                }
            }
        }
    }

   
    virtual ~ScaffoldedSupertree(){}
    ScaffoldedSupertree()
        :TaxonomyDependentTreeProcessor<TreeMappedWithSplits>(),
         numErrors(0),
         doReportAllContested(false),
         doConstructSupertree(false),
         taxonomyAsSource(nullptr),
         currDotFileIndex(0),
         debuggingOutput(false), 
         emitScaffoldDotFiles(false) {
    }

    void reportAllConflicting(std::ostream & out, bool verbose) {
        std::map<std::size_t, unsigned long> nodeMappingDegree;
        std::map<std::size_t, unsigned long> passThroughDegree;
        std::map<std::size_t, unsigned long> loopDegree;
        unsigned long totalContested = 0;
        unsigned long redundContested = 0;
        unsigned long totalNumNodes = 0;
        for (auto nd : iter_node_internal(*taxonomy)) {
            const auto & thr = _getEmdeddingForNode(nd);
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
                const auto & thr = _getEmdeddingForNode(nd);
                std::vector<NodeWithSplits *> aliasedBy = getNodesAliasedBy(nd, *taxonomy);
                thr.reportIfContested(out, nd, treePtrByIndex, aliasedBy, otCLI.verbose);
            }
        }
        for (auto tr : idListForDotExport) {
            auto nd = taxonomy->getData().getNodeForOttId(tr);
            if (nd == nullptr) {
                throw OTCError(std::string("Unrecognized OTT ID in list of OTT IDs to export to DOT: ") + std::to_string(tr));
            }
            for (auto n : iter_pre_n_const(nd)) {
                const auto & thr = _getEmdeddingForNode(n);
                writeDOTExport(out, thr, n, treePtrByIndex, false, false);
            }
        }
        return true;
    }

    bool processTaxonomyTree(OTCLI & otCLI) override {
        debuggingOutput = otCLI.verbose;
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::processTaxonomyTree(otCLI);
        checkTreeInvariants(*taxonomy);
        suppressMonotypicTaxaPreserveDeepestDangle(*taxonomy);
        checkTreeInvariants(*taxonomy);
        for (NodeWithSplits * nd : iter_node(*taxonomy)) {
            _getEmdeddingForNode(nd);
        }
        otCLI.getParsingRules().setOttIdForInternals = false;
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
        embedNewTree(*taxonomy, *raw, treeIndex);
        otCLI.err << "# pathPairings = " << pathPairings.size() << '\n';
        return true;
    }

    bool cloneTaxonomyAsASourceTree() {
        assert(taxonomy != nullptr);
        assert(taxonomyAsSource == nullptr);
        std::unique_ptr<TreeMappedWithSplits> tree = std::move(cloneTree(*taxonomy));
        taxonomyAsSource = tree.get();
        std::size_t treeIndex = inputTreesToIndex.size();
        inputTreesToIndex[std::move(tree)] = treeIndex;
        treePtrByIndex.push_back(taxonomyAsSource);
        // suppress the internal node OTT IDs from the des
        OttIdSet internalIDs;
        for (auto nd : iter_post_internal(*taxonomyAsSource)) {
            if (nd->hasOttId()) {
                internalIDs.insert(nd->getOttId());
            }
            auto & d = nd->getData().desIds;
            for (auto o : internalIDs) {
                d.erase(o);
            }
        }
        // Store the tree's filename
        taxonomyAsSource->setName("TAXONOMY");
        embedTaxonomyClone(*taxonomy, *taxonomyAsSource, treeIndex);
        return true;
    }
};

bool handleReportAllFlag(OTCLI & otCLI, const std::string &);
bool handleReportOnNodesFlag(OTCLI & otCLI, const std::string &);
bool handleDotNodesFlag(OTCLI & otCLI, const std::string &narg);
bool handleSuperTreeFlag(OTCLI & otCLI, const std::string &narg);

bool handleReportAllFlag(OTCLI & otCLI, const std::string &) {
    ScaffoldedSupertree * proc = static_cast<ScaffoldedSupertree *>(otCLI.blob);
    assert(proc != nullptr);
    proc->doReportAllContested = true;
    return true;
}
bool handleSuperTreeFlag(OTCLI & otCLI, const std::string &) {
    ScaffoldedSupertree * proc = static_cast<ScaffoldedSupertree *>(otCLI.blob);
    assert(proc != nullptr);
    proc->doConstructSupertree = true;
    return true;
}

bool handleReportOnNodesFlag(OTCLI & otCLI, const std::string &narg) {
    ScaffoldedSupertree * proc = static_cast<ScaffoldedSupertree *>(otCLI.blob);
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
bool handleOttForestDOTFlag(OTCLI & otCLI, const std::string &narg) {
    long conv = -1;
    if (!char_ptr_to_long(narg.c_str(), &conv) || conv < 0) {
        throw OTCError(std::string("Expecting a positive number as an ott ID after -z flag. Offending word: ") + narg);
    }
    ottIDBeingDebugged = conv;
    return true;
}
bool handleOttScaffoldDOTFlag(OTCLI & otCLI, const std::string &) {
    ScaffoldedSupertree * proc = static_cast<ScaffoldedSupertree *>(otCLI.blob);
    assert(proc != nullptr);
    proc->emitScaffoldDotFiles = true;
    return true;
}
bool handleDotNodesFlag(OTCLI & otCLI, const std::string &narg) {
    ScaffoldedSupertree * proc = static_cast<ScaffoldedSupertree *>(otCLI.blob);
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
    ScaffoldedSupertree proc;
    otCLI.addFlag('a',
                  "Write a report of all contested nodes",
                  handleReportAllFlag,
                  false);
    otCLI.addFlag('s',
                  "Compute a supertree",
                  handleSuperTreeFlag,
                  false);
    otCLI.addFlag('b',
                  "ARG should be a list of OTT IDs. A status report will be generated for those nodes",
                  handleReportOnNodesFlag,
                  true);
    otCLI.addFlag('d',
                  "ARG should be a list of OTT IDs. A DOT file of the nodes will be generated ",
                  handleDotNodesFlag,
                  true);
    otCLI.addFlag('z',
                  "ARG should be an OTT ID. A series DOT files will be generated for the forest created during the resolution of this OTT ID ",
                  handleOttForestDOTFlag,
                  true);
    otCLI.addFlag('y',
                  "requests DOT export of the embedded tree during the supertree operation - only for use on small examples!",
                  handleOttScaffoldDOTFlag,
                  false);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}

