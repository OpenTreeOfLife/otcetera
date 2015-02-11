#include "otc/newick.h"
#include "otc/test_harness.h"
using namespace otc;


char testSingleCharLabelPoly(const TestHarness &);
char testWordLabelPoly(const TestHarness &);
char testUnderscores(const TestHarness &);
char testWhitespaceHandling(const TestHarness &);
char testQuotedWordLabelPoly(const TestHarness &);
char testCommentPoly(const TestHarness &);
char testUnbalanced(const TestHarness &);
char testUnbalancedToManyClose(const TestHarness &);


char genericTokenTest(const TestHarness &th, const std::string &fn, const std::vector<std::string> & expected);
char genericOTCParsingErrorTest(const TestHarness &th, const std::string &fn);

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

char genericOTCParsingErrorTest(const TestHarness &th, const std::string &fn) {
	std::ifstream inp;
	if (!th.openTestFile(fn, inp)) {
		return 'U';
	}
	NewickTokenizer tokenizer(inp, th.getFilePath(fn));
	try {
		for (auto token : tokenizer) {
		}
	} catch (const OTCParsingError &) {
		return '.';
	}
	std::cerr << "No OTCParsingError was raised!\n";
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

char testQuotedWordLabelPoly(const TestHarness &th) {
	const std::vector<std::string> expected = {"(", "AB", ",", "BC", ",", "CD", ")", ";"};
	return genericTokenTest(th, "quotedwords-polytomy.tre", expected);
}

char testCommentPoly(const TestHarness &th) {
	const std::vector<std::string> expected = {"(", "AB", ",", "BC", ",", "CD", ")", ";"};
	return genericTokenTest(th, "polytomy-with-comments.tre", expected);
}
char testUnderscores(const TestHarness &th) {
	const std::vector<std::string> expected = {"(", "A B", ",", "B C", ",", "C_D", ")", ";"};
	return genericTokenTest(th, "underscore-handling.tre", expected);
}
char testWhitespaceHandling(const TestHarness &th) {
	const std::vector<std::string> expected = {"(", "A B",
		                                       ",", "B C",
		                                       ",", "C_D",
		                                       ",", "EX FZ", //warn but tolerate unquoted..
		                                       ",", "G HY", //warn but normalize unquoted..
		                                       ",", "I J", //warn but normalize unquoted even with newlines, I guess...
		                                       ")", ";"};
	return genericTokenTest(th, "whitespace-handling.tre", expected);
}
char testUnbalanced(const TestHarness &th) {
	return genericOTCParsingErrorTest(th, "unbalanced.tre");
}
char testUnbalancedToManyClose(const TestHarness &th) {
	return genericOTCParsingErrorTest(th, "unbalanced-too-many-closed.tre");
}

int main(int argc, char *argv[]) {
	TestHarness th(argc, argv);
	TestsVec tests{TestFn("testSingleCharLabelPoly", testSingleCharLabelPoly)
				   , TestFn("testWordLabelPoly", testWordLabelPoly)
				   , TestFn("testQuotedWordLabelPoly", testQuotedWordLabelPoly)
				   , TestFn("testCommentPoly", testCommentPoly)
				   , TestFn("testUnderscores", testUnderscores)
				   , TestFn("testUnbalanced", testUnbalanced)
				   , TestFn("testUnbalancedToManyClose", testUnbalancedToManyClose)
				   , TestFn("testWhitespaceHandling", testWhitespaceHandling)
				   
				  };
	return th.runTests(tests);
}

