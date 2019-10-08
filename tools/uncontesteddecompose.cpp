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
        assert(!scaffoldNd->is_tip());
        auto & thr = _get_embedding_for_node(scaffoldNd);

        LOG(INFO) << " exportOrCollapse for ott" << scaffoldNd->get_ott_id() << " outdegree = " << scaffoldNd->get_out_degree() << " numLoopTrees = " << thr.get_num_loop_trees() << " numLoops = " << thr.get_total_num_loops();
        if (thr.is_contested()) {
            //auto p = scaffoldNd->get_parent();
            LOG(INFO) << "    Contested";
            if (documentP != nullptr) {
                auto treeIndToContestingNodeMap = thr.get_how_tree_contests_monophyly_maps();
                json treeIDToNodeMapJSON;
                for (auto& [treei, parToChildSetMap] : treeIndToContestingNodeMap)
                {
                    json parToChildSetJSON = json::array();
                    for (auto& [parNode, childSet] : parToChildSetMap)
                    {
                        json childSetAsJSONList = json::array();
                        for (auto childP : childSet)
                        {
                            childSetAsJSONList.push_back(childP->get_name());
                        }

                        json pcsObj = {{"parent",parNode->get_name()},{"children_from_taxon", childSetAsJSONList}};
                        parToChildSetJSON.push_back(pcsObj);
                    }

                    const auto ct = treePtrByIndex.at(treei);
                    treeIDToNodeMapJSON[ct->get_name()] = parToChildSetJSON;
                }
                std::string ottIdStr = "ott" + std::to_string(scaffoldNd->get_ott_id());
                (*documentP)[ottIdStr] = treeIDToNodeMapJSON;
            }
            thr.collapse_group(*scaffoldNd, sc);
        } else {
            //thr.debug_node_embeddings(" focal node before export", false, scaffoldNdToNodeEmbedding);
            //if (scaffoldNd->get_parent()) {
            //    _get_embedding_for_node(scaffoldNd->get_parent()).debug_node_embeddings(" parent before export", true, scaffoldNdToNodeEmbedding);
            //}
            LOG(INFO) << "    Uncontested";
            auto fn = thr.export_subproblem_and_resolve(*scaffoldNd, exportDir, exportStream, sc);
            //if (scaffoldNd->get_parent()) {
            //    _get_embedding_for_node(scaffoldNd->get_parent()).debug_node_embeddings("after export", true, scaffoldNdToNodeEmbedding);
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
            sc.prune_tips_mapped_to_contested_taxa = false;
        }
        std::list<NodeWithSplits * > postOrder;
        for (auto nd : iter_post(*taxonomy)) {
            if (nd->is_tip()) {
                assert(nd->has_ott_id());
                // this is only needed for monotypic cases in which a tip node
                //  may have multiple OTT Ids in its des_ids set
                _get_embedding_for_node(nd).set_ott_id_for_exit_embeddings(nd,
                                                                   nd->get_ott_id(),
                                                                   scaffoldNdToNodeEmbedding);
            } else {
                postOrder.push_back(nd);
            }
            //_get_embedding_for_node(nd).debug_node_embeddings(" getting postorder", true, scaffoldNdToNodeEmbedding);
        }
        for (auto nd : postOrder) {
            assert(!nd->is_tip());
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
        clone_taxonomy_as_a_source_tree();
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
    otCLI.add_flag('e',
                  "ARG should be the name of a directory. A .tre file will be written to that directory for each subproblem",
                  handleExportSubproblems,
                  true);
    otCLI.add_flag('o',
                  "If present, the trees will be exported to standard output",
                  handleExportToStdoutSubproblems,
                  false);
    otCLI.add_flag('r',
                  "If present, the tips in input trees which are mapped to contested taxa. The default behavior is to prune these tips",
                  handleRetainTipsMapToContestedTaxaSubproblems,
                  false);
    otCLI.add_flag('x',
                  "ARG should be a file path. a line listing the name (but not the full path) over every created .tre file will be written to this file.",
                  handleListSubproblemIds,
                  true);
    otCLI.add_flag('c',
                  "ARG should be a file path. A JSON representation of the trees that contest each taxon will be written to that filepath.",
                  handleContestingLog,
                  true);
    return tax_dependent_tree_processing_main(otCLI, argc, argv, proc, 2, true);
}

