#include "otc/embedding_cli.h"
using namespace otc;

class UncontestedTaxonDecompose : public EmbeddingCLI {
    public:
    std::string exportDir;
    std::string subproblemIdFile;
    std::ostream * exportStream;
    bool userRequestsRetentionOfTipsMappedToContestedTaxa;

    virtual ~UncontestedTaxonDecompose(){}
    UncontestedTaxonDecompose()
        :EmbeddingCLI(),
        exportStream(nullptr),
        userRequestsRetentionOfTipsMappedToContestedTaxa(false) {
    }

    void exportOrCollapse(NodeWithSplits * scaffoldNd, SupertreeContextWithSplits & sc) {
        assert(!scaffoldNd->isTip());
        auto & thr = _getEmbeddingForNode(scaffoldNd);

        LOG(INFO) << " exportOrCollapse for ott" << scaffoldNd->getOttId() << " outdegree = " << scaffoldNd->getOutDegree() << " numLoopTrees = " << thr.getNumLoopTrees() << " numLoops = " << thr.getTotalNumLoops();
        if (thr.isContested()) {
            //thr.debugNodeEmbedding("before thr.constructPhyloGraphAndCollapseIfNecessary", true, scaffoldNdToNodeEmbedding);
            //auto p = scaffoldNd->getParent();
            LOG(INFO) << "    Contested";
            thr.constructPhyloGraphAndCollapseIfNecessary(*scaffoldNd, sc);
            //_getEmbeddingForNode(p).debugNodeEmbedding("after thr.constructPhyloGraphAndCollapseIfNecessary", true, scaffoldNdToNodeEmbedding);
        } else {
            //thr.debugNodeEmbedding(" focal node before export", false, scaffoldNdToNodeEmbedding);
            //if (scaffoldNd->getParent()) {
            //    _getEmbeddingForNode(scaffoldNd->getParent()).debugNodeEmbedding(" parent before export", true, scaffoldNdToNodeEmbedding);
            //}
            LOG(INFO) << "    Uncontested";
            thr.exportSubproblemAndResolve(*scaffoldNd, exportDir, exportStream, sc);
            //if (scaffoldNd->getParent()) {
            //    _getEmbeddingForNode(scaffoldNd->getParent()).debugNodeEmbedding("after export", true, scaffoldNdToNodeEmbedding);
            //}
        }
    }

    void exportSubproblems(OTCLI &) {
        TreeMappedWithSplits * tax = taxonomy.get();
        SupertreeContextWithSplits sc{treePtrByIndex, scaffoldNdToNodeEmbedding, *tax};
        if (userRequestsRetentionOfTipsMappedToContestedTaxa) {
            sc.pruneTipsMappedToContestedTaxa = false;
        }
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
            //_getEmbeddingForNode(nd).debugNodeEmbedding(" getting postorder", true, scaffoldNdToNodeEmbedding);
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
bool handleExportToStdoutSubproblems(OTCLI & otCLI, const std::string &narg);
bool handleRetainTipsMapToContestedTaxaSubproblems(OTCLI & otCLI, const std::string &narg);
bool handleListSubproblemIds(OTCLI & otCLI, const std::string &narg);

bool handleExportToStdoutSubproblems(OTCLI & otCLI, const std::string &) {
    UncontestedTaxonDecompose * proc = static_cast<UncontestedTaxonDecompose *>(otCLI.blob);
    assert(proc != nullptr);
    proc->exportStream = &(otCLI.out);
    return true;
}

bool handleRetainTipsMapToContestedTaxaSubproblems(OTCLI & otCLI, const std::string &) {
    UncontestedTaxonDecompose * proc = static_cast<UncontestedTaxonDecompose *>(otCLI.blob);
    assert(proc != nullptr);
    proc->userRequestsRetentionOfTipsMappedToContestedTaxa = true;
    return true;
}

bool handleExportSubproblems(OTCLI & otCLI, const std::string &narg) {
    UncontestedTaxonDecompose * proc = static_cast<UncontestedTaxonDecompose *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a directory after the -e argument.");
    }
    proc->exportDir = narg;
    return true;
}

bool handleListSubproblemIds(OTCLI & otCLI, const std::string &narg) {
    UncontestedTaxonDecompose * proc = static_cast<UncontestedTaxonDecompose *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a filepath after the -x argument.");
    }
    proc->subproblemIdFile = narg;
    return true;
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-uncontested-decompose",
                "takes at least 2 newick file paths: a full taxonomy tree, and some number of input trees, and -e flag to specify an export directory",
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
    otCLI.addFlag('r',
                  "If present, the tips in input trees which are mapped to contested taxa. The default behavior is to prune these tips",
                  handleRetainTipsMapToContestedTaxaSubproblems,
                  false);
    otCLI.addFlag('x',
                  "ARG should be the name of a file. a line listing the name (but not the full path) over every created .tre file will be written to this file.",
                  handleListSubproblemIds,
                  true);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}

