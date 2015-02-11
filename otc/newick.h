#ifndef OTCETERA_NEWICK_H
#define OTCETERA_NEWICK_H
#include <iostream>
#include <fstream>
#include <stack>
#include <stdexcept>
#include "otc/otc_base_includes.h"
#include "otc/tree.h"
#include "otc/error.h"

namespace otc {

typedef std::shared_ptr<const std::string> ConstStrPtr;
struct FilePosStruct {
	FilePosStruct()
		:pos(0U),
		lineNumber(0U),
		colNumber(0U),
		filepath(ConstStrPtr(nullptr)) {
	}
	FilePosStruct(const ConstStrPtr fp)
		:pos(0U),
		lineNumber(0U),
		colNumber(0U),
		filepath(fp) {
	}
	FilePosStruct(std::size_t position, std::size_t line, std::size_t column, const ConstStrPtr fp)
		:pos(position),
		lineNumber(line),
		colNumber(column),
		filepath(fp) {
	}
	std::string describe() const {
		std::string message;
		message = "At line ";
		message += std::to_string(1 + this->lineNumber);
		message += ", column ";
		message += std::to_string(1 + this->colNumber);
		message += ", filepos ";
		message += std::to_string(1 + this->pos);
		try {
			if (filepath) {
				message += " of \"";
				message += *filepath;
				message += "\"";
			} else {
				message += " of <UNKNOWN FILENAME>";
			}
		} catch (...) {
		}
		return message;
	}

	std::size_t pos;
	std::size_t lineNumber;
	std::size_t colNumber;
	ConstStrPtr filepath;
};

class OTCParsingError: public OTCError {
	public:
		OTCParsingError(const char * msg, char offending, const FilePosStruct &position)
			:OTCError(msg),
			frag(msg),
			offendingChar(offending),
			pos(position) {
				message = this->generate_message();
			}
		std::string generate_message() const noexcept {
			try {
				std::string m = "Error found \"";
				if (this->offendingChar != '\0') {
					m += this->offendingChar;
				}
				m += "\" ";
				m += this->frag;
				m += " ";
				m += this->pos.describe();
				return m;
			} catch (...) {// to guarantee noexcept...
			}
			return frag;
		}

	private:
		const std::string frag;
		const char offendingChar;
		const FilePosStruct pos;
};
class NewickTokenizer {
	public:
		enum newick_token_state_t {
				NWK_NOT_IN_TREE, // before first open-parens
				NWK_OPEN, // token was open-parens
				NWK_CLOSE, // token was close-parens
				NWK_COMMA, // token was ,
				NWK_COLON, // token was :
				NWK_BRANCH_INFO, // token was text after :
				NWK_LABEL, // token was text (but not after a :)
				NWK_SEMICOLON
			};
		NewickTokenizer(std::istream &inp, const std::string & filepath)
			:inputStream(inp),
			inpFilepath(nullptr) {
			if (filepath.length() > 0) {
				const std::string * raw = new std::string(filepath);
				inpFilepath = std::shared_ptr<const std::string>(raw);
			}
		}
		class iterator;
		class Token {
			public:
				const std::string & content() const {
					return this->tokenContent;
				}
				const std::vector<std::string> & commentVec() const {
					return this->comments;
				}
			private:
				Token(const std::string &content,
					  const FilePosStruct & startPosition,
					  const FilePosStruct & endPosition,
					  const std::vector<std::string> &embeddedComments,
					  newick_token_state_t tokenState)
					:tokenContent(content), 
					startPos(startPosition),
					endPos(endPosition),
					comments(embeddedComments),
					state(tokenState) {
					LOG(TRACE) << "created token for \"" << content << "\"";
				}
			public:
				const std::string tokenContent;
				const FilePosStruct startPos;
				const FilePosStruct endPos;
				const std::vector<std::string> comments;
				const newick_token_state_t state;
				
				friend class NewickTokenizer::iterator;
		};

		class iterator : std::forward_iterator_tag {
			public:
				bool operator==(const iterator & other) const {
					//LOG(TRACE) << "Equality test atEnd = " << this->atEnd << " other.atEnd = " << other.atEnd;
					if (other.atEnd) {
						return this->atEnd;
					}
					if (this->atEnd) {
						return false;
					}
					return (this->currentPos == other.currentPos) && (&(this->inputStream) == &(other.inputStream));
				}
				bool operator!=(const iterator & other) const {
					//LOG(TRACE) << "Inequality test atEnd = " << this->atEnd << " other.atEnd = " << other.atEnd << '\n';
					return !(*this == other);
				}
				Token operator*() const {
					LOG(TRACE) << "* operator";
					return Token(currWord, *prevPos, *currentPos, comments, currTokenState);
				}
				iterator & operator++() {
					LOG(TRACE) << "increment";
					if (this->atEnd) {
						throw std::out_of_range("Incremented a dead NewickTokenizer::iterator");
					}
					if (this->currTokenState == NWK_SEMICOLON) {
						this->atEnd = true;
						this->resetToken();
					} else {
						this->resetToken();
						this->consumeNextToken();
					}
					return *this;
				}
			private:
				void consumeNextToken();
				bool advanceToNextNonWhitespace(char &);
				void finishReadingComment();
				void finishReadingUnquoted(bool continuingLabel);
				void finishReadingQuotedStr();
				void onLabelExit(char nextChar, bool enteringFromWhitespace);
				char peek() {
					if (!pushed.empty()) {
						return pushed.top();
					}
					char c = (char) this->inputStream.rdbuf()->sgetc();
					pushed.push(c);
					return c;
				}
				void push(char c) {
					this->pushed.push(c);
					if (c == '\n') {
						this->currentPos->colNumber = lastLineInd;
						this->currentPos->lineNumber -= 1;
					} else {
						this->currentPos->colNumber -= 1;
					}
					this->currentPos->pos -= 1;
				}
				void throwSCCErr(char c) const __attribute__ ((noreturn));
				//deals with \r\n as \n Hence "LogicalChar"
				bool advanceReaderOneLogicalChar(char & c) {
					if (!pushed.empty()) {
						c = pushed.top();
						pushed.pop();
					} else {
						if (this->atEnd) {
							c = EOF;
							return false;
						}
						c = (char) (this->inputStream.rdbuf())->sbumpc();
					}
					if (c == EOF) {
						this->atEnd = true;
					} else {
						this->currentPos->pos += 1;
						if (13 == c || 10 == c) {
							if (13 == c ) { // deal with \r\n as a newline
								if (this->inputStream.rdbuf()->sgetc() == 10) {//peeks at the next char
									(inputStream.rdbuf())->sbumpc();
									this->currentPos->pos += 1;
								}
							}
							c = '\n';
							this->lastLineInd = this->currentPos->colNumber;
							this->currentPos->colNumber = 0;
							this->currentPos->lineNumber += 1;
						} else {
							this->currentPos->colNumber += 1;
						}
					}
					return !this->atEnd;
				}
				void resetToken() {
					this->currWord.clear();
					this->prevTokenState = this->currTokenState;
					std::swap(currentPos, prevPos);
					currentPos->pos = prevPos->pos;
					currentPos->lineNumber = prevPos->lineNumber;
					currentPos->colNumber = prevPos->colNumber;
					comments.clear();
				}
				iterator(std::istream &inp, ConstStrPtr filepath)
					:inputStream(inp),
					inpFilepath(filepath),
					atEnd(!inp.good()),
					firstPosSlot(filepath),
					secondPosSlot(filepath),
					currTokenState(NWK_NOT_IN_TREE),
					prevTokenState(NWK_NOT_IN_TREE),
					numUnclosedParens(0) {
					currentPos = &firstPosSlot;
					prevPos = &secondPosSlot;
					LOG(TRACE) << "create live";
					++(*this);
				}
				iterator(std::istream &inp) // USE in end() ONLY!
					:inputStream(inp),
					inpFilepath(nullptr),
					atEnd(true),
					firstPosSlot(nullptr),
					secondPosSlot(nullptr),
					currTokenState(NWK_NOT_IN_TREE),
					prevTokenState(NWK_NOT_IN_TREE),
					numUnclosedParens(0) {
					LOG(TRACE) << "create dead";
					
				}
				std::istream & inputStream;
				ConstStrPtr inpFilepath;
				bool atEnd;
				FilePosStruct firstPosSlot;
				FilePosStruct secondPosSlot;
				FilePosStruct * currentPos; // alias
				FilePosStruct * prevPos;    // alias
				std::string currWord;
				newick_token_state_t currTokenState;
				newick_token_state_t prevTokenState;
				long numUnclosedParens;
				std::stack<char> pushed;
				std::vector<std::string> comments;
				std::size_t lastLineInd;
				friend class NewickTokenizer;
		};
		iterator begin() {
			iterator b(this->inputStream, this->inpFilepath);
			return b;
		}
		iterator end() {
			LOG(TRACE) << "NewickTokenizer.end()";
			return iterator(this->inputStream);
		}

	private:
		std::istream & inputStream;
		ConstStrPtr inpFilepath;
};

//Takes wide istream and (optional) filepath (just used for error reporting if not empty)
template<typename T, typename U>
std::unique_ptr<RootedTree<T, U> > readNextNewick(std::istream &inp, const std::string & filepath);


void newickParseNodeInfo(RootedTree<RTNodeNoData, RTreeNoData> & ,
						 RootedTreeNode<RTNodeNoData> &,
						 const NewickTokenizer::Token * labelToken,
						 const NewickTokenizer::Token * colonToken, // used for comment
						 const NewickTokenizer::Token * branchLenToken);
void newickCloseNodeHook(RootedTree<RTNodeNoData, RTreeNoData> & ,
						 RootedTreeNode<RTNodeNoData> &,
						 const NewickTokenizer::Token & closeToken);

inline void newickCloseNodeHook(RootedTree<RTNodeNoData, RTreeNoData> & ,
								RootedTreeNode<RTNodeNoData> & ,
								const NewickTokenizer::Token & ) {
}

inline void newickParseNodeInfo(RootedTree<RTNodeNoData, RTreeNoData> & ,
								RootedTreeNode<RTNodeNoData> & node,
								const NewickTokenizer::Token * labelToken,
								const NewickTokenizer::Token * , // used for comment
								const NewickTokenizer::Token * ) {
	if (labelToken) {
		node.SetName(labelToken->content());
	}
}

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
