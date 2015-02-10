#ifndef OTCETERA_NEWICK_H
#define OTCETERA_NEWICK_H
#include <iostream>
#include <fstream>
#include "otc/otcetera.h"
#include "otc/tree.h"

namespace otc {

/* TREE IO */
#if 0
// non wide string versions
template<typename T>
std::unique_ptr<RootedTree<T> > readNextNewick(std::istream &inp);
template<typename T>
unsigned int readNewickStream(std::istream &inp, bool (*callback)(std::unique_ptr<RootedTree<T> >));

template<typename T>
inline std::unique_ptr<RootedTree<T> > readNextNewick(std::istream &inp) {
	return std::unique_ptr<RootedTree<T> > (new RootedTree<T>());
}

template<typename T>
inline unsigned int readNewickStream(std::istream &inp, bool (*callback)(std::unique_ptr<RootedTree<T> >)) {
	auto c = 0U;
	for (;;) {
		std::unique_ptr<RootedTree<T> > nt = readNextNewick<T>(inp);
		if (nt == nullptr) {
			return c;
		}
		++c;
		auto cbr = callback(std::move(nt));
		if (!cbr) {
			return c;
		}
	}
}
#endif

typedef std::shared_ptr<const std::string> ConstStrPtr;
struct FilePosStruct {
	FilePosStruct()
		:pos(0U),
		lineNumber(0U),
		colNumber(0U),
		filepath(nullptr) {
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
		message += 1 + this->lineNumber;
		message += ", column ";
		message += 1 + this->colNumber;
		message += ", filepos ";
		message += 1 + this->pos;
		if (filepath) {
			message += " of \"";
			message += filepath
			message += "\"";
		} else {
			message += " of <UNKNOWN FILENAME>";
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
			:message(msg),
			offendingChar(offending),
			pos(position) {
			}
		const char * what () const noexcept {
			try {
				m = "Error found \"";
				m += offendingChar;
				m += "\" ";
				m += this.message;
				m += this->pos.describe();
				return m.c_str();
			} catch (...) { // to guarantee noexcept...
			}
			return message.c_str();
		}

	private:
		const std::string message;
		const char offendingChar;
		const FilePosStruct pos;
		const m;
};
class NewickTokenizer {
	public:
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
			private:
				Token(const std::string &content,
					  const FilePosStruct & startPosition,
					  const FilePosStruct & endPosition)
					:tokenContent(content), 
					startPos(startPosition),
					endPos(endPosition) {
				}
				const std::string tokenContent;
				const FilePosStruct startPos;
				const FilePosStruct endPos;
				friend class NewickTokenizer::iterator;
		};

		class iterator : std::forward_iterator_tag {
			enum newick_token_state_t {
				NWK_NOT_IN_TREE, // before first open-parens
				NWK_OPEN, // token was open-parens
				NWK_CLOSE, // token was close-parens
				NWK_COMMA, // token was ,
				NWK_COLON, // token was :
				NWK_BRANCH_INFO, // token was text after :
				NWK_LABEL // token was text (but not after a :) 
			};
			public:
				bool operator==(const iterator & other) const {
					LOG(TRACE) << "Equality test atEnd = " << this->atEnd << " other.atEnd = " << other.atEnd << '\n';
					if (other.atEnd) {
						return this->atEnd;
					}
					if (this->atEnd) {
						return false;
					}
					return (this->currentPos == other.currentPos) && (&(this->inputStream) == &(other.inputStream));
				}
				bool operator!=(const iterator & other) const {
						LOG(TRACE) << "Inequality test atEnd = " << this->atEnd << " other.atEnd = " << other.atEnd << '\n';
						return !(*this == other);
					}
				Token operator*() const {
					LOG(TRACE) << "* operator\n";
					return Token(currWord, *prevPos, *currentPos);
				}
				iterator & operator++() {
					LOG(TRACE) << "increment\n";
					if (this->atEnd) {
						throw std::out_of_range("Incremented a dead NewickTokenizer::iterator");
					}
					this->resetToken();
					this->consumeNextToken();
					return *this;
				}
			private:
				void consumeNextToken();
				//deals with \r\n as \n Hence "LogicalChar"
				bool advanceReaderOneLogicalChar(char & c) {
					c = (char) (this->inputStream.rdbuf())->sbumpc();
					if (c == EOF) {
						this->atEnd = true;
					} else {
						this->currentPos->pos += 1;
						if (13 == c || 10 == c) {
							if (13 == c ) { // deal with \r\n as a newline
								if (this->inputStream.rdbuf())->sgetc() == 10) {//peeks at the next char
									(inputStream.rdbuf())->sbumpc();
									this->currentPos->pos += 1;
								}
							}
							c = '\n';
							this->currentPos->colNumber = 0;
							this->currentPos->lineNumber += 1;
						} else {
							this->currentPos->colNumber += 1;
						}
					}
				}
				void resetToken() {
					this->currWord.clear();
					this->prevTokenState = this->currTokenState;
					std::swap(currentPos, prevPos);
					currentPos->pos = prevPos->pos;
					currentPos->lineNumber = prevPos->lineNumber;
					currentPos->colNumber = prevPos->colNumber;
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
					LOG(TRACE) << "create live\n";
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
					LOG(TRACE) << "create dead\n";
					
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
				friend class NewickTokenizer;
		};
		iterator begin() {
			iterator b(this->inputStream, this->inpFilepath);
			return b;
		}
		iterator end() {
			return iterator(this->inputStream);
		}

	private:
		std::istream & inputStream;
		ConstStrPtr inpFilepath;
};

//Takes wide istream and (optional) filepath (just used for error reporting if not empty)
template<typename T>
std::unique_ptr<RootedTree<T> > readNextNewick(std::istream &inp, const std::string & filepath);

template<typename T>
inline std::unique_ptr<RootedTree<T> > readNextNewick(std::istream &inp, const std::string & filepath) {
	assert(inp.good());
	NewickTokenizer tokenizer(inp, filepath);
	auto i = 0U;
	for (auto token : tokenizer) {
		std::cout << "token = \"" << token.content() << "\"\n"; 
		if (i > 1000U) {
			break;
		}
	}
	return std::unique_ptr<RootedTree<T> > (new RootedTree<T>());
}



}// namespace otc
#endif
