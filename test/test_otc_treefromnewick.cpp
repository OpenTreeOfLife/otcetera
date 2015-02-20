#include "otc/newick.h"
#include "otc/util.h"
#include "otc/test_harness.h"
using namespace otc;

class TestValidTreeStruct {
		const std::string filename;
	public:
		TestValidTreeStruct(const std::string & fn)
			:filename(fn) {
		}
		char runTest(const TestHarness &h) const {
			auto fp = h.getFilePath(filename);
			//std::cerr << "fn = " << fp<< '\n';
			std::ifstream inp;
			if (!openUTF8File(fp, inp)) {
				return 'U';
			}
			for (;;) {
				auto nt = readNextNewick<RTNodeNoData, RTreeNoData>(inp, filename);
				if (nt == nullptr) {
					return '.';
				}
			}
		}
};


int main(int argc, char *argv[]) {
	std::vector<std::string> validfilenames = {"abc-newick.tre", 
						   "words-polytomy.tre", 
						   "bifurcating.tre",
						   "monotypic.tre", 
						   "branch-lengths.tre",
						   "quotedwords-polytomy.tre",
						   "polytomy-with-comments.tre",
						   "underscore-handling.tre", 
						   "whitespace-handling.tre"};
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

