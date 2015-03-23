#include "otc/embedding_cli.h"
using namespace otc;

class UncontestedTaxonDecompose : public EmbeddingCLI {
    public:
    std::string exportDir;
    std::ostream * exportStream;

    virtual ~UncontestedTaxonDecompose(){}
    UncontestedTaxonDecompose()
        :EmbeddingCLI(),
        exportStream(nullptr) {
    }

    void exportOrCollapse(NodeWithSplits * scaffoldNd, SupertreeContextWithSplits & sc) {
        assert(!scaffoldNd->isTip());
        auto & thr = _getEmbeddingForNode(scaffoldNd);

        LOG(INFO) << " exportOrCollapse for ott" << scaffoldNd->getOttId() << " outdegree = " << scaffoldNd->getOutDegree() << " numLoopTrees = " << thr.getNumLoopTrees() << " numLoops = " << thr.getTotalNumLoops();
        if (thr.isContested()) {
            thr.debugNodeEmbedding(true, scaffoldNdToNodeEmbedding);
            LOG(INFO) << "    Contested";
            thr.constructPhyloGraphAndCollapseIfNecessary(*scaffoldNd, sc);
            _getEmbeddingForNode(scaffoldNd->getParent()).debugNodeEmbedding(true, scaffoldNdToNodeEmbedding);
        } else {
            thr.debugNodeEmbedding(false, scaffoldNdToNodeEmbedding);
            if (scaffoldNd->getParent()) {
                _getEmbeddingForNode(scaffoldNd->getParent()).debugNodeEmbedding(true, scaffoldNdToNodeEmbedding);
            }
            LOG(INFO) << "    Uncontested";
            thr.exportSubproblemAndFakeResolution(*scaffoldNd, exportDir, exportStream, sc);
            if (scaffoldNd->getParent()) {
                _getEmbeddingForNode(scaffoldNd->getParent()).debugNodeEmbedding(true, scaffoldNdToNodeEmbedding);
            }
        }
    }

    void exportSubproblems(OTCLI &) {
        TreeMappedWithSplits * tax = taxonomy.get();
        SupertreeContextWithSplits sc{treePtrByIndex, scaffoldNdToNodeEmbedding, *tax};
        std::list<NodeWithSplits * > postOrder;
        for (auto nd : iter_post(*taxonomy)) {
            if (nd->isTip()) {
                assert(nd->hasOttId());
                // this is only needed for monotypic cases in which a tip node
                //  may have multiple OTT Ids in its desIds set
                _getEmbeddingForNode(nd).setOttIdForExitEmbeddings(nd,
                                                                   nd->getOttId(),
                                                                   scaffoldNdToNodeEmbedding);
            } else {
                postOrder.push_back(nd);
            }
            _getEmbeddingForNode(nd).debugNodeEmbedding(true, scaffoldNdToNodeEmbedding);
        }
        for (auto nd : postOrder) {
            assert(!nd->isTip());
            exportOrCollapse(nd, sc);
            for (auto nnd : postOrder) {
                _getEmbeddingForNode(nnd).debugNodeEmbedding(true, scaffoldNdToNodeEmbedding);
            }
        }
    }

    bool summarize(OTCLI &otCLI) override {
        cloneTaxonomyAsASourceTree();
        exportSubproblems(otCLI);
        return true;
    }
};

bool handleExportSubproblems(OTCLI & otCLI, const std::string &narg);
bool handleExportToStdoutSubproblems(OTCLI & otCLI, const std::string &narg);

bool handleExportToStdoutSubproblems(OTCLI & otCLI, const std::string &) {
    UncontestedTaxonDecompose * proc = static_cast<UncontestedTaxonDecompose *>(otCLI.blob);
    assert(proc != nullptr);
    proc->exportStream = &(otCLI.out);
    return true;
}

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
    otCLI.addFlag('o',
                  "If present, the trees will be exported to standard output",
                  handleExportToStdoutSubproblems,
                  false);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}

