#include "otc/embedding_cli.h"
#include "json.hpp"
using namespace otc;
using json = nlohmann::json;

class UncontestedTaxonDecompose : public EmbeddingCLI {
    public:
    std::string exportDir;
    std::string subproblemIdFile;
    std::string contestingLogFile;
    std::ostream * exportStream;
    std::ostream * subproblemIdStream;
    bool userRequestsRetentionOfTipsMappedToContestedTaxa;

    virtual ~UncontestedTaxonDecompose(){}
    UncontestedTaxonDecompose()
        :EmbeddingCLI(),
        exportStream(nullptr),
        subproblemIdStream(nullptr),
        userRequestsRetentionOfTipsMappedToContestedTaxa(false) {
    }

    void exportOrCollapse(NodeWithSplits * scaffoldNd, SupertreeContextWithSplits & sc, json * documentP) {
        assert(!scaffoldNd->isTip());
        auto & thr = _get_embedding_for_node(scaffoldNd);

        LOG(INFO) << " exportOrCollapse for ott" << scaffoldNd->get_ott_id() << " outdegree = " << scaffoldNd->getOutDegree() << " numLoopTrees = " << thr.getNumLoopTrees() << " numLoops = " << thr.getTotalNumLoops();
        if (thr.isContested()) {
            //auto p = scaffoldNd->getParent();
            LOG(INFO) << "    Contested";
            if (documentP != nullptr) {
                auto treeIndToContestingNodeMap = thr.getHowTreeContestsMonophylyMaps();
                json treeIDToNodeMapJSON;
                for (auto treeCNMPair : treeIndToContestingNodeMap) {
                    auto & treei = treeCNMPair.first;
                    auto & parToChildSetMap = treeCNMPair.second;
                    const auto ct = treePtrByIndex.at(treei);
                    json parToChildSetJSON = json::array();
                    for (auto parChildSetPair : parToChildSetMap) {
                        json pcsObj;
                        const auto & parNode = parChildSetPair.first;
                        const auto & childSet = parChildSetPair.second;
                        json childSetAsJSONList = json::array();
                        for (auto childP : childSet) {
                            childSetAsJSONList.push_back(childP->get_name());
                        }
                        pcsObj["parent"] = parNode->get_name();
                        pcsObj["children_from_taxon"] = childSetAsJSONList;
                        parToChildSetJSON.push_back(pcsObj);
                    }
                    treeIDToNodeMapJSON[ct->get_name()] = parToChildSetJSON;
                }
                std::string ottIdStr = "ott" + std::to_string(scaffoldNd->get_ott_id());
                (*documentP)[ottIdStr] = treeIDToNodeMapJSON;
            }
            thr.collapseGroup(*scaffoldNd, sc);
        } else {
            //thr.debugNodeEmbedding(" focal node before export", false, scaffoldNdToNodeEmbedding);
            //if (scaffoldNd->getParent()) {
            //    _get_embedding_for_node(scaffoldNd->getParent()).debugNodeEmbedding(" parent before export", true, scaffoldNdToNodeEmbedding);
            //}
            LOG(INFO) << "    Uncontested";
            auto fn = thr.exportSubproblemAndResolve(*scaffoldNd, exportDir, exportStream, sc);
            //if (scaffoldNd->getParent()) {
            //    _get_embedding_for_node(scaffoldNd->getParent()).debugNodeEmbedding("after export", true, scaffoldNdToNodeEmbedding);
            //}
            if ((subproblemIdStream != nullptr) && (!fn.empty())) {
                *subproblemIdStream << fn << '\n';
            }
        }
    }

    void exportSubproblems(OTCLI &, json * documentP) {
        TreeMappedWithSplits * tax = taxonomy.get();
        SupertreeContextWithSplits sc{treePtrByIndex, scaffoldNdToNodeEmbedding, *tax};
        if (userRequestsRetentionOfTipsMappedToContestedTaxa) {
            sc.pruneTipsMappedToContestedTaxa = false;
        }
        std::list<NodeWithSplits * > postOrder;
        for (auto nd : iter_post(*taxonomy)) {
            if (nd->isTip()) {
                assert(nd->has_ott_id());
                // this is only needed for monotypic cases in which a tip node
                //  may have multiple OTT Ids in its desIds set
                _get_embedding_for_node(nd).set_ott_idForExitEmbeddings(nd,
                                                                   nd->get_ott_id(),
                                                                   scaffoldNdToNodeEmbedding);
            } else {
                postOrder.push_back(nd);
            }
            //_get_embedding_for_node(nd).debugNodeEmbedding(" getting postorder", true, scaffoldNdToNodeEmbedding);
        }
        for (auto nd : postOrder) {
            assert(!nd->isTip());
            exportOrCollapse(nd, sc, documentP);
        }
    }

    bool summarize(OTCLI &otCLI) override {
        std::ofstream sif;
        if (!subproblemIdFile.empty()) {
            sif.open(subproblemIdFile.c_str());
            if (!sif.good()) {
                throw OTCError("Could not open subproblem ID file");
            }
            subproblemIdStream = &sif;
        }
        json * documentP = nullptr;
        json document;
        std::ofstream clf;
        if (!contestingLogFile.empty()) {
            clf.open(contestingLogFile.c_str());
            if (!clf.good()) {
                throw OTCError("Could not open contesting log file");
            }
            documentP = &document;
        }
        cloneTaxonomyAsASourceTree();
        exportSubproblems(otCLI, documentP);
        subproblemIdStream = nullptr;
        if (documentP != nullptr) {
            clf << document.dump(1) << std::endl;
        }
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

bool handleContestingLog(OTCLI & otCLI, const std::string &narg) {
    UncontestedTaxonDecompose * proc = static_cast<UncontestedTaxonDecompose *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a filepath after the -c argument.");
    }
    proc->contestingLogFile = narg;
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
                  "ARG should be a file path. a line listing the name (but not the full path) over every created .tre file will be written to this file.",
                  handleListSubproblemIds,
                  true);
    otCLI.addFlag('c',
                  "ARG should be a file path. A JSON representation of the trees that contest each taxon will be written to that filepath.",
                  handleContestingLog,
                  true);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}

