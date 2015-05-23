#include "otc/newick.h"
#include "otc/util.h"
#include "otc/test_harness.h"
using namespace otc;

typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;
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
            ConstStrPtr filenamePtr = ConstStrPtr(new std::string(filename));
            FilePosStruct pos(filenamePtr);
            for (;;) {
                ParsingRules pr;
                auto nt = readNextNewick<Tree_t>(inp, pos, pr);
                return (nt != nullptr ? '.': 'F');
            }
        }
};

int main(int argc, char *argv[]) {
    std::vector<std::string> validfilenames = {"noids-abcnewick.tre", 
                           "noids-wordspolytomy.tre", 
                           "noids-bifurcating.tre",
                           "noids-monotypic.tre", 
                           "noids-branchlengths.tre",
                           "noids-quotedwordspolytomy.tre",
                           "noids-polytomywithcomments.tre",
                           "noids-underscorehandling.tre", 
                           "noids-whitespacehandling.tre"};
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

