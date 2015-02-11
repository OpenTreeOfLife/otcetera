#include "otc/newick.h"
#include "otc/test_harness.h"
using namespace otc;


char testSingleCharLabelPoly(const TestHarness &);
char testSingleCharLabelPoly(const TestHarness &th) {
	const std::string fn = "abc-newick.tre";
	std::ifstream inp;
	if (!th.openTestFile(fn, inp)) {
		return 'U';
	}
	const std::vector<std::string> expected = {"(", "A", ",", "B", ",", "C", ")", ";"};
	std::vector<std::string> obtained;
	NewickTokenizer tokenizer(inp, th.getFilePath(fn));
	for (auto token : tokenizer) {
		obtained.push_back(token.content());
		}
	if (testVecElementEquality(expected, obtained)) {
		return '.';
	}
	return 'F';
}

int main(int argc, char *argv[]) {
	TestHarness th(argc, argv);
	TestsVec tests{TestFn("testSingleCharLabelPoly", testSingleCharLabelPoly)
				  };
	return th.runTests(tests);
}

