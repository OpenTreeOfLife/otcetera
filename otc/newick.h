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
template<typename T>
std::unique_ptr<T> readNextNewick(std::istream &inp, FilePosStruct & pos, const ParsingRules &parsingRules);

template<typename T>
inline std::unique_ptr<T> readNextNewick(std::istream &inp, FilePosStruct & pos, const ParsingRules &parsingRules) {
    assert(inp.good());
    NewickTokenizer tokenizer(inp, pos);
    auto tokenIt = tokenizer.begin();
    if (tokenIt == tokenizer.end()) {
        return std::unique_ptr<T>(nullptr);
    }
    std::stack<typename T::node_type *> nodeStack;
    T * rawTreePtr = new T();
    std::unique_ptr<T> treePtr(rawTreePtr);
    typename T::node_type * currNode = rawTreePtr->createRoot();
    // If we read a label or colon, we might consume multiple tokens;
    for (; tokenIt != tokenizer.end(); ) {
        const NewickTokenizer::Token topOfLoopToken = *tokenIt;
        if (topOfLoopToken.state == NewickTokenizer::NWK_OPEN) {
            nodeStack.push(currNode);
            currNode = rawTreePtr->createChild(currNode);
            ++tokenIt;
        } else if (topOfLoopToken.state == NewickTokenizer::NWK_CLOSE) {
            assert(!nodeStack.empty()); // NewickTokenizer should thrown an exception if unbalanced
            newickCloseNodeHook(*rawTreePtr, *currNode, topOfLoopToken, parsingRules);
            currNode = nodeStack.top();
            nodeStack.pop();
            ++tokenIt;
        } else if (topOfLoopToken.state == NewickTokenizer::NWK_COMMA) {
            assert(!nodeStack.empty()); // NewickTokenizer should thrown an exception if unbalanced
            newickCloseNodeHook(*rawTreePtr, *currNode, topOfLoopToken, parsingRules);
            currNode = rawTreePtr->createSib(currNode);
            ++tokenIt;
        } else if (topOfLoopToken.state == NewickTokenizer::NWK_LABEL) {
            ++tokenIt;
            const NewickTokenizer::Token colonToken = *tokenIt;
            if (colonToken.state == NewickTokenizer::NWK_COLON) {
                ++tokenIt;
                const NewickTokenizer::Token brLenToken = *tokenIt;
                assert(brLenToken.state == NewickTokenizer::NWK_BRANCH_INFO);
                newickParseNodeInfo(*rawTreePtr, *currNode, &topOfLoopToken, &colonToken, &brLenToken, parsingRules);
                ++tokenIt;
            } else {
                newickParseNodeInfo(*rawTreePtr, *currNode, &topOfLoopToken, nullptr, nullptr, parsingRules);
            }
        } else if (topOfLoopToken.state == NewickTokenizer::NWK_COLON) {
            ++tokenIt;
            const NewickTokenizer::Token brLenToken = *tokenIt;
            assert(brLenToken.state == NewickTokenizer::NWK_BRANCH_INFO);
            newickParseNodeInfo(*rawTreePtr, *currNode, nullptr, &topOfLoopToken, &brLenToken, parsingRules);
            ++tokenIt;
        } else {
            assert(topOfLoopToken.state == NewickTokenizer::NWK_SEMICOLON);
            break;
        }
    }
    postParseHook(*treePtr, parsingRules);
    pos.setLocationInFile(tokenIt.getCurrPos());
    return treePtr;
}

}// namespace otc
#endif
