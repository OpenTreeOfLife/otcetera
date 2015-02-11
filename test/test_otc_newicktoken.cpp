#include "otc/newick.h"
#include "otc/test_harness.h"
using namespace otc;


char testSingleCharLabelPoly(const TestHarness &);
char testWordLabelPoly(const TestHarness &);

char genericTokenTest(const TestHarness &th, const std::string &fn, const std::vector<std::string> & expected);

char genericTokenTest(const TestHarness &th, const std::string &fn, const std::vector<std::string> & expected){
	std::ifstream inp;
	if (!th.openTestFile(fn, inp)) {
		return 'U';
	}
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
char testSingleCharLabelPoly(const TestHarness &th) {
	const std::vector<std::string> expected = {"(", "A", ",", "B", ",", "C", ")", ";"};
	return genericTokenTest(th, "abc-newick.tre", expected);
}

char testWordLabelPoly(const TestHarness &th) {
	const std::vector<std::string> expected = {"(", "AB", ",", "BC", ",", "CD", ")", ";"};
	return genericTokenTest(th, "words-polytomy.tre", expected);
}

int main(int argc, char *argv[]) {
	TestHarness th(argc, argv);
	TestsVec tests{TestFn("testSingleCharLabelPoly", testSingleCharLabelPoly),
				   TestFn("testWordLabelPoly", testWordLabelPoly)
				  };
	return th.runTests(tests);
}

