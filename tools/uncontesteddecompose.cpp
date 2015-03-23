#include "otc/embeddingCLI.h"
using namespace otc;

class UncontestedTaxonDecompose : public EmbeddingCLI {
    public:
    std::string exportDir;

    virtual ~UncontestedTaxonDecompose(){}
    UncontestedTaxonDecompose()
        :EmbeddingCLI() {
    }

    void exportOrCollapse(NodeWithSplits * scaffoldNd, SupertreeContextWithSplits & sc) {
        assert(!scaffoldNd->isTip());
        auto & thr = _getEmdeddingForNode(scaffoldNd);
        LOG(INFO) << "  outdegree = " << scaffoldNd->getOutDegree() << " numLoopTrees = " << thr.getNumLoopTrees() << " numLoops = " << thr.getTotalNumLoops();
        if (thr.isContested()) {
            LOG(INFO) << "    Contested";
            if (thr.highRankingTreesPreserveMonophyly(sc.numTrees)) {
                thr.resolveGivenContestedMonophyly(*scaffoldNd, sc);
            } else {
                thr.constructPhyloGraphAndCollapseIfNecessary(*scaffoldNd, sc);
            }
        } else {
            LOG(INFO) << "    Uncontested";
            if (exportDir.empty()) {
                thr.resolveGivenUncontestedMonophyly(*scaffoldNd, sc);
            } else {
                thr.exportSubproblemAndFakeResolution(*scaffoldNd, exportDir, sc);
            }
        }
    }

    void exportSubproblems(OTCLI &) {
        TreeMappedWithSplits * tax = taxonomy.get();
        SupertreeContextWithSplits sc{treePtrByIndex, scaffoldNdToNodeEmbedding, *tax};
        std::list<NodeWithSplits * > postOrder;
        for (auto nd : iter_post_internal(*taxonomy)) {
            assert(!nd->isTip());
            postOrder.push_back(nd);
        }
        for (auto nd : postOrder) {
            assert(!nd->isTip());
            exportOrCollapse(nd, sc);
        }
    }

    bool summarize(OTCLI &otCLI) override {
        cloneTaxonomyAsASourceTree();
        exportSubproblems(otCLI);
        return true;
    }
};

bool handleExportSubproblems(OTCLI & otCLI, const std::string &narg);

bool handleExportSubproblems(OTCLI & otCLI, const std::string &narg) {
    UncontestedTaxonDecompose * proc = static_cast<UncontestedTaxonDecompose *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a list of IDs after the -b argument.");
    }
    proc->exportDir = narg;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otcuncontesteddecompose",
                "takes at least 2 newick file paths: a full taxonomy tree, and some number of input trees, and -e flag to specify an export directory.",
                "taxonomy.tre inp1.tre inp2.tre");
    UncontestedTaxonDecompose proc;
    otCLI.addFlag('e',
                  "ARG should be the name of a directory. A .tre file will be written to that directory for each subproblem",
                  handleExportSubproblems,
                  true);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}

