#include "otc/newick.h"
#include "otc/util.h"
#include "otc/test_harness.h"
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
            std::cerr << "tv.size() == " << tv.size() << "\n";
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

