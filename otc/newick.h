#ifndef OTCETERA_NEWICK_H
#define OTCETERA_NEWICK_H
#include <iostream>
#include <fstream>
#include <stack>
#include <stdexcept>
#include <sstream>
#include "otc/otc_base_includes.h"
#include "otc/tree.h"
#include "otc/newick_tokenizer.h"
#include "otc/parse_newick_data.h"
#include "otc/error.h"

namespace otc {
//Takes wide istream and (optional) filepath (just used for error reporting if not empty)
template<typename T>
std::unique_ptr<T> read_next_newick(std::istream &inp, FilePosStruct & pos, const ParsingRules &parsingRules);

template<typename T>
inline std::unique_ptr<T> read_next_newick(std::istream &inp, FilePosStruct & pos, const ParsingRules &parsingRules) {
    assert(inp.good());
    NewickTokenizer tokenizer(inp, pos);
    auto tokenIt = tokenizer.begin();
    if (tokenIt == tokenizer.end()) {
        return std::unique_ptr<T>(nullptr);
    }
    std::stack<typename T::node_type *> nodeStack;
    T * rawTreePtr = new T();
    std::unique_ptr<T> treePtr(rawTreePtr);
    typename T::node_type * currNode = rawTreePtr->create_root();
    // If we read a label or colon, we might consume multiple tokens;
    for (; tokenIt != tokenizer.end(); ) {
        const NewickTokenizer::Token topOfLoopToken = *tokenIt;
        if (topOfLoopToken.state == NewickTokenizer::NWK_OPEN) {
            nodeStack.push(currNode);
            currNode = rawTreePtr->create_child(currNode);
            ++tokenIt;
        } else if (topOfLoopToken.state == NewickTokenizer::NWK_CLOSE) {
            assert(!nodeStack.empty()); // NewickTokenizer should thrown an exception if unbalanced
            newick_close_node_hook(*rawTreePtr, *currNode, topOfLoopToken, parsingRules);
            currNode = nodeStack.top();
            nodeStack.pop();
            ++tokenIt;
        } else if (topOfLoopToken.state == NewickTokenizer::NWK_COMMA) {
            assert(!nodeStack.empty()); // NewickTokenizer should thrown an exception if unbalanced
            newick_close_node_hook(*rawTreePtr, *currNode, topOfLoopToken, parsingRules);
            currNode = rawTreePtr->create_sib(currNode);
            ++tokenIt;
        } else if (topOfLoopToken.state == NewickTokenizer::NWK_LABEL) {
            ++tokenIt;
            const NewickTokenizer::Token colonToken = *tokenIt;
            if (colonToken.state == NewickTokenizer::NWK_COLON) {
                ++tokenIt;
                const NewickTokenizer::Token brLenToken = *tokenIt;
                assert(brLenToken.state == NewickTokenizer::NWK_BRANCH_INFO);
                newick_parse_node_info(*rawTreePtr, *currNode, &topOfLoopToken, &colonToken, &brLenToken, parsingRules);
                ++tokenIt;
            } else {
                newick_parse_node_info(*rawTreePtr, *currNode, &topOfLoopToken, nullptr, nullptr, parsingRules);
            }
        } else if (topOfLoopToken.state == NewickTokenizer::NWK_COLON) {
            ++tokenIt;
            const NewickTokenizer::Token brLenToken = *tokenIt;
            assert(brLenToken.state == NewickTokenizer::NWK_BRANCH_INFO);
            newick_parse_node_info(*rawTreePtr, *currNode, nullptr, &topOfLoopToken, &brLenToken, parsingRules);
            ++tokenIt;
        } else {
            assert(topOfLoopToken.state == NewickTokenizer::NWK_SEMICOLON);
            break;
        }
    }
    post_parse_hook(*treePtr, parsingRules);
    pos.set_location_in_file(tokenIt.get_curr_pos());
    return treePtr;
}

template<typename T>
inline std::unique_ptr<T> tree_from_newick_string(const std::string& s, const ParsingRules & Rules) {
    std::istringstream sfile(s);
    FilePosStruct pos;
    auto tree = read_next_newick<T>(sfile, pos, Rules);
    if (not tree)
	throw OTCParsingError("Newick string is empty (no tokens)");
    return tree;
}

template<typename T>
inline std::unique_ptr<T> tree_from_newick_string(const std::string& s) {
    ParsingRules rules;
    rules.require_ott_ids = false;
    return tree_from_newick_string<T>(s,rules);
}

template <typename T>
inline std::unique_ptr<T> first_newick_tree_from_file(const std::string& filename, const ParsingRules& Rules)
{
    std::ifstream inp;
    if (!open_utf8_file(filename, inp)) {
	throw OTCError()<<"Could not open \""<<filename<<"\"";
    }
    LOG(INFO) << "reading \"" << filename << "\"...";
    ConstStrPtr filenamePtr = ConstStrPtr(new std::string(filename));
    FilePosStruct pos(filenamePtr);
    return read_next_newick<T>(inp, pos, Rules);
}

template <typename T>
inline std::unique_ptr<T> first_newick_tree_from_file(const std::string& filename)
{
    ParsingRules rules;
    rules.require_ott_ids = false;
    return first_newick_tree_from_file<T>(filename, rules);
}

}// namespace otc


#endif
