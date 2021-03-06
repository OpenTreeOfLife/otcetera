#include "otc/newick.h"
#include "otc/util.h"
#include "otc/test_harness.h"
#include <sstream>
using namespace otc;

template<typename T>
inline std::string getNewick(const T *nd) {
    std::ostringstream out;
    write_newick(out, nd);
    out << ";\n";
    return out.str();
}

typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;
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
            for (;;) {
                ParsingRules pr;
                auto nt = read_next_newick<Tree_t>(inp, pos, pr);
                if (filename == "3genus-synth.tre") {
                    const std::string expected[] = {
                        "(((A1_ott1,(A2_ott2,A3_ott3))A_ott4,B1_ott5),B2_ott6,B3_ott7,((C1_ott9,C2_ott10),C3_ott11)C_ott12)life_ott14;\n",
                        "((A1_ott1,(A2_ott2,A3_ott3))A_ott4,B1_ott5);\n",
                        "(A1_ott1,(A2_ott2,A3_ott3))A_ott4;\n",
                        "A1_ott1;\n",
                        "(A2_ott2,A3_ott3);\n",
                        "A2_ott2;\n",
                        "A3_ott3;\n",
                        "B1_ott5;\n",
                        "B2_ott6;\n",
                        "B3_ott7;\n",
                        "((C1_ott9,C2_ott10),C3_ott11)C_ott12;\n",
                        "(C1_ott9,C2_ott10);\n",
                        "C1_ott9;\n",
                        "C2_ott10;\n",
                        "C3_ott11;\n"
                    };
                    return doIterTest(expected, *nt);
                }
                if (filename == "3genus-taxonomy.tre") {
                    const std::string expected[] = {
                        "((A1_ott1,A2_ott2,A3_ott3)A_ott4,(B1_ott5,B2_ott6,B3_ott7)B_ott8,(C1_ott9,C2_ott10,C3_ott11)C_ott12)life_ott14;\n",
                        "(A1_ott1,A2_ott2,A3_ott3)A_ott4;\n",
                        "A1_ott1;\n",
                        "A2_ott2;\n",
                        "A3_ott3;\n",
                        "(B1_ott5,B2_ott6,B3_ott7)B_ott8;\n",
                        "B1_ott5;\n",
                        "B2_ott6;\n",
                        "B3_ott7;\n",
                        "(C1_ott9,C2_ott10,C3_ott11)C_ott12;\n",
                        "C1_ott9;\n",
                        "C2_ott10;\n",
                        "C3_ott11;\n"
                    };
                    return doIterTest(expected, *nt);
                }
            }
        }
        char doIterTest(const std::string expected [], Tree_t &tree) const {
            auto i = 0;
            for (auto nd: iter_pre_n_const(tree.get_root())) {
                if (expected[i] != getNewick(nd)) {
                    std::cerr << expected[i] << " != " << getNewick(nd) << '\n';
                    return 'F';
                } else {
                }
                ++i;
            }
            return '.';
        }
};

int main(int argc, char *argv[]) {
    std::vector<std::string> validfilenames = {"3genus-synth.tre", "3genus-taxonomy.tre"};
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

