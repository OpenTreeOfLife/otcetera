#include "otc/newick.h"
#include "otc/util.h"
#include "otc/test_harness.h"
#include "otc/tree_data.h"
#include "otc/node_embedding.h"
#include "otc/tree_iter.h"
#include "otc/greedy_forest.h"
#include "otc/embedded_tree.h"
using namespace otc;

typedef TreeMappedWithSplits Tree_t;
class TestValidTreeStruct {
        const std::string filename;
    public:
        TestValidTreeStruct(const std::string & fn)
            :filename(fn) {
        }
        char runTest(const TestHarness &h) const {
            auto fp = h.get_filepath(filename);
            std::ifstream inp;
            if (!open_utf8_file(fp, inp)) {
                return 'U';
            }
            ConstStrPtr filenamePtr = ConstStrPtr(new std::string(filename));
            FilePosStruct pos(filenamePtr);
            std::list<std::unique_ptr<TreeMappedWithSplits> > tv;
            for (;;) {
                ParsingRules pr;
                auto nt = read_next_newick<Tree_t>(inp, pos, pr);
                if (nt == nullptr) {
                    break;
                }
                tv.push_back(std::move(nt));
            }
            EmbeddedTree et;
            Tree_t fakeScaffold;
            auto r = fakeScaffold.create_root();
            auto & fsd = fakeScaffold.get_data().ott_id_to_node;
            std::set<long> ids;
            for (const auto & tp : tv) {
                const Tree_t & tree = *tp;
                const OttIdSet & td = tree.get_root()->get_data().des_ids;
                for (auto oid : td) {
                    if (!contains(ids, oid)) {
                        fsd[oid] = fakeScaffold.create_child(r);
                        ids.insert(oid);
                    }
                }
            }
            Tree_t fakePhylo;
            std::vector<Tree_t *> ftv;
            ftv.push_back(&fakePhylo);

            SupertreeContextWithSplits sc{ftv,
                                          et._get_scaffold_nd_to_node_embedding(),
                                          fakeScaffold};
            GreedyBandedForest<NodeWithSplits, NodeWithSplits> gpf{-1};
            int treeInd = 0;
            long groupInd = 0;
            for (const auto & tp : tv) {
                const Tree_t & tree = *tp;
                const OttIdSet * incGroup = nullptr;
                for (auto nd : iter_child(*tree.get_root())) {
                    if (!nd->is_tip()) {
                        if (incGroup != nullptr) {
                            return 'F';
                        }
                        incGroup = &(nd->get_data().des_ids);
                    }
                }
                if (incGroup == nullptr) {
                    return 'F';
                }
                const OttIdSet & leaf_set = tree.get_root()->get_data().des_ids;
                std::cerr << groupInd + 1 << '\n';
                gpf.attempt_to_add_grouping(*incGroup, leaf_set, treeInd, groupInd++, nullptr);
            }
            NodeEmbeddingWithSplits emptyEmbedding(r);
            gpf.finish_resolution_of_embedded_clade(*r, &emptyEmbedding, &sc);
            return '.';
        }
};

int main(int argc, char *argv[]) {
    std::vector<std::string> validfilenames = { //"multi-band-simple-phylo-statements.tre" 
                                                //};
                                                /*
                                                , "multi-band-phylo-statements.tre" 
                                                , "shortest-crash-phylo-statements-small-tax.tre"
                                               ,"shortest-crash-phylo-statements.tre"
                                               , "shorter-phylo-statements.tre"
                                               , */
                                               //"phylo-statements.tre"
                                               };
    TestHarness th(argc, argv);
    TestsVec tests;
    for (auto fn : validfilenames) {
        //const TestValidTreeStruct tvts(fn);
        const TestValidTreeStruct tvts{fn};
        TestCallBack tcb = [tvts](const TestHarness &h) {
            return tvts.runTest(h);
        };
        const TestFn tf{fn, tcb};
        tests.push_back(tf);
    }
    return th.run_tests(tests);
}

