#include "otc/newick.h"
#include "otc/util.h"
#include "otc/test_harness.h"
#include "otc/tree_data.h"
#include "otc/embedding.h"
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
            auto fp = h.getFilePath(filename);
            std::ifstream inp;
            if (!openUTF8File(fp, inp)) {
                return 'U';
            }
            std::list<std::unique_ptr<TreeMappedWithSplits> > tv;
            for (;;) {
                ParsingRules pr;
                auto nt = readNextNewick<Tree_t>(inp, filename, pr);
                if (nt == nullptr) {
                    break;
                }
                tv.push_back(std::move(nt));
            }
            if (tv.size() != 164) {
              return 'F';
            }
            EmbeddedTree et;
            Tree_t fakeScaffold;
            auto r = fakeScaffold.createRoot();
            auto & fsd = fakeScaffold.getData().ottIdToNode;
            std::set<long> ids;
            for (const auto & tp : tv) {
                const Tree_t & tree = *tp;
                const OttIdSet & td = tree.getRoot()->getData().desIds;
                for (auto oid : td) {
                    if (!contains(ids, oid)) {
                        fsd[oid] = fakeScaffold.createChild(r);
                        ids.insert(oid);
                    }
                }
            }

            SupertreeContextWithSplits sc{1, et.taxoToEmbedding, fakeScaffold};
            GreedyPhylogeneticForest<NodeWithSplits, NodeWithSplits> gpf{-1};
            int treeInd = 0;
            long groupInd = 0;
            for (const auto & tp : tv) {
                const Tree_t & tree = *tp;
                const OttIdSet * incGroup = nullptr;
                for (auto nd : iter_child(*tree.getRoot())) {
                    if (!nd->isTip()) {
                        if (incGroup != nullptr) {
                            return 'F';
                        }
                        incGroup = &(nd->getData().desIds);
                    }
                }
                if (incGroup == nullptr) {
                    return 'F';
                }
                const OttIdSet & leafSet = tree.getRoot()->getData().desIds;
                gpf.attemptToAddGrouping(nullptr, *incGroup, leafSet, treeInd, groupInd++, nullptr);
            }
            NodeEmbeddingWithSplits emptyEmbedding(r);
            gpf.finishResolutionOfEmbeddedClade(*r, &emptyEmbedding, &sc);
            return '.';
        }
};

int main(int argc, char *argv[]) {
    std::vector<std::string> validfilenames = {"phylo-statements.tre"};
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
    return th.runTests(tests);
}

