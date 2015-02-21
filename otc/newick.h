#ifndef OTCETERA_NEWICK_H
#define OTCETERA_NEWICK_H
#include <iostream>
#include <fstream>
#include <stack>
#include <stdexcept>
#include "otc/otc_base_includes.h"
#include "otc/tree.h"
#include "otc/newick_tokenizer.h"
#include "otc/parse_newick_data.h"
#include "otc/error.h"

namespace otc {


//Takes wide istream and (optional) filepath (just used for error reporting if not empty)
template<typename T, typename U>
std::unique_ptr<RootedTree<T, U> > readNextNewick(std::istream &inp, const std::string & filepath);

template<typename T, typename U>
inline std::unique_ptr<RootedTree<T, U> > readNextNewick(std::istream &inp, const std::string & filepath) {
	assert(inp.good());
	NewickTokenizer tokenizer(inp, filepath);
	auto tokenIt = tokenizer.begin();
	if (tokenIt == tokenizer.end()) {
		return std::unique_ptr<RootedTree<T, U> >(nullptr);
	}
	std::stack<RootedTreeNode<T> *> nodeStack;
	RootedTree<T, U> * rawTreePtr = new RootedTree<T, U>();
	std::unique_ptr<RootedTree<T, U> > treePtr(rawTreePtr);
	RootedTreeNode<T> * currNode = rawTreePtr->CreateRoot();
	// If we read a label or colon, we might consume multiple tokens;
	for (;tokenIt != tokenizer.end(); ) {
		const NewickTokenizer::Token topOfLoopToken = *tokenIt;
		if (topOfLoopToken.state == NewickTokenizer::NWK_OPEN) {
			nodeStack.push(currNode);
			currNode = rawTreePtr->CreateChild(currNode);
			++tokenIt;
		} else if (topOfLoopToken.state == NewickTokenizer::NWK_CLOSE) {
			assert(!nodeStack.empty()); // NewickTokenizer should thrown an exception if unbalanced
			newickCloseNodeHook(*rawTreePtr, *currNode, topOfLoopToken);
			currNode = nodeStack.top();
			nodeStack.pop();
			++tokenIt;
		} else if (topOfLoopToken.state == NewickTokenizer::NWK_COMMA) {
			assert(!nodeStack.empty()); // NewickTokenizer should thrown an exception if unbalanced
			currNode = rawTreePtr->CreateSib(currNode);
			++tokenIt;
		} else if (topOfLoopToken.state == NewickTokenizer::NWK_LABEL) {
			++tokenIt;
			const NewickTokenizer::Token colonToken = *tokenIt;
			if (colonToken.state == NewickTokenizer::NWK_COLON) {
				++tokenIt;
				const NewickTokenizer::Token brLenToken = *tokenIt;
				assert(brLenToken.state == NewickTokenizer::NWK_BRANCH_INFO);
				newickParseNodeInfo(*rawTreePtr, *currNode, &topOfLoopToken, &colonToken, &brLenToken);
				++tokenIt;
			} else {
				newickParseNodeInfo(*rawTreePtr, *currNode, &topOfLoopToken, nullptr, nullptr);
			}
		} else if (topOfLoopToken.state == NewickTokenizer::NWK_COLON) {
			++tokenIt;
			const NewickTokenizer::Token brLenToken = *tokenIt;
			assert(brLenToken.state == NewickTokenizer::NWK_BRANCH_INFO);
			newickParseNodeInfo(*rawTreePtr, *currNode, nullptr, &topOfLoopToken, &brLenToken);
		} else {
			assert(topOfLoopToken.state == NewickTokenizer::NWK_SEMICOLON);
			break;
		}
	}
	return treePtr;
}



}// namespace otc
#endif
