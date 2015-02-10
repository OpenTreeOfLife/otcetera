#include "otc/newick.h"
#include <cstring>
#include <string>
namespace otc {

static const char * _EARLY_SEMICOLON = "Unexpected ; with open parentheses not balanced";
static const char * _ILL_AFTER_CLOSE = "Illegal character after \")\" character. Expecting \",\" or a label or a colon.";
static const char * _ILL_AFTER_LABEL = "Illegal character after label. Expecting ( or a label";
static const char * _ILL_AFTER_BRANCH_INFO = "Illegal character after branch info. Expecting ( or a label";
static const char * _ILL_AFTER_OPEN = "Illegal character after \"(\" character. Expecting ( or a label";
static const char * _ILL_AFTER_COMMA = "Illegal character after \",\" character. Expecting ( or a label.";
static const char * _ILL_AFTER_COLON = "Illegal character after \":\" character. Expecting a branch length.";
static const char * _ILL_NO_SEMICOLON = "Expecting ; after a newick description.";
static const char * _ILL_FIRST_CHAR = "Expecting a newick tree to start with \"(\"";


void NewickTokenizer::iterator::onLabelExit(char n) {
	bool whitespaceFound = false;
	if (std::strchr("(),:;", n) == nullptr) {
		if (!std::isgraph(n)) {
			whitespaceFound = true;
			if (!advanceToNextNonWhitespace(n)) {
				return;
			}
			if (std::strchr("(),:;", n) != nullptr) { // exit stripping whitespace
				this->push(n);
				return;
			}
			this->currWord += ' '; // replace any series of one white space with one " ". Not sure what to do, this seems reasonable.
		}
		if (n == '[') {
			this->advanceReaderOneLogicalChar(n);
			assert(n == '[');
			finishReadingComment();
			n = this->peek();
			onLabelExit(n);
			return;
		}
		LOG(WARNING) << "Unexpected continuation of a label after a quoted string in newick. Next character is: \"" << n << "\"";
		if (n == '\'') {
			this->advanceReaderOneLogicalChar(n);
			assert(n == '\'');
			finishReadingQuotedStr();
			return;
		}
		this->advanceReaderOneLogicalChar(n);
		this->finishReadingUnquoted(true);
		return;
	}
}
void NewickTokenizer::iterator::finishReadingQuotedStr(){
	for (;;) {
		char c;
		if(!advanceReaderOneLogicalChar(c)) {
			throw OTCParsingError("Unexpected EOF in quoted string", '\0', (*this->currentPos));
		}
		if (c == '\'') {
			auto n = this->peek();
			if (n == '\'') {
				this->currWord += c;
			} else {
				this->onLabelExit(n);
				return;
			}
		} else {
			this->currWord += c;
		}
	}
}
void NewickTokenizer::iterator::finishReadingUnquoted(bool continuingLabel){
	for (;;) {
		char c;
		if(!advanceReaderOneLogicalChar(c)) {
			throw OTCParsingError("Unexpected EOF in label. Expecting a ; to end a newick", '\0', (*this->currentPos));
		}
		if (!isgraph(c)) {
			if (!advanceToNextNonWhitespace(c)) {
				throw OTCParsingError("Unexpected EOF in label. Expecting a ; to end a newick", '\0', (*this->currentPos));
			}
			if (continuingLabel) {
				LOG(WARNING) << "Whitespace found in unquoted label - will be converted to a single space";
			}
			this->push(c);
			this->onLabelExit(c);
			return;
		}
		if (c == '\'') {
			LOG(WARNING) << "single-quoted string found in unquoted label";
			this->finishReadingQuotedStr();
			return;
		} else if (c == '[') {
			finishReadingComment();
		} else if (c == '_') {
			this->currWord += ' ';
		} else {
			this->currWord += c;
		}
	}
}
void NewickTokenizer::iterator::finishReadingComment(){
	auto numOpenComments = 1U;
	std::string comment;
	for (;;) {
		char c;
		if(!advanceReaderOneLogicalChar(c)) {
			throw OTCParsingError("Unexpected EOF in comment", '\0', (*this->currentPos));
		}
		if (c == ']') {
			numOpenComments -= 1;
			if (numOpenComments == 0) {
				this->comments.push_back(comment);
			}
			return;
		}
		if (c == '[') {
			numOpenComments += 1;
		}
		comment += c;
	}

}
bool NewickTokenizer::iterator::advanceToNextNonWhitespace(char & c){
	for (;;) {
		if(!advanceReaderOneLogicalChar(c)) {
			return false;
		}
		if (std::isgraph(c)) {
			return true;
		}
	}
}

void NewickTokenizer::iterator::throwSCCErr(char n) const {
	if (this->prevTokenState == NWK_OPEN) {
		throw OTCParsingError(_ILL_AFTER_OPEN, n, *this->currentPos);
	}
	if (this->prevTokenState == NWK_COLON) {
		throw OTCParsingError(_ILL_AFTER_COLON, n, *this->currentPos);
	}
	if (this->prevTokenState == NWK_COMMA) {
		throw OTCParsingError(_ILL_AFTER_OPEN, n, *this->currentPos);
	}
	assert(false);
	throw OTCError("should be unreachable");
}

void NewickTokenizer::iterator::consumeNextToken() {
	char n;
	// a loop to accumulate comments
	while (true) {
		if (advanceToNextNonWhitespace(n)) {
			if (this->prevTokenState == NWK_NOT_IN_TREE) {
				if (n == '(') {
					this->currTokenState = NWK_OPEN;
					this->numUnclosedParens += 1;
					this->currWord.assign(1, '(');
					return;
				} else if (n == '[') {
					this->finishReadingComment();
				} else {
					throw OTCParsingError(_ILL_FIRST_CHAR, n, *this->currentPos);
				}
			}
			else {
				if (std::strchr("(),:;[\'", n) == nullptr) {
					this->currTokenState = (this->prevTokenState == NWK_COLON ? NWK_BRANCH_INFO : NWK_LABEL);
					this->currWord.assign(1, n);
					this->finishReadingUnquoted(false);
					return;
				} else {
					if (n != ';' && this->numUnclosedParens <= 0) {
						assert(this->numUnclosedParens == 0);
						throw OTCParsingError(_ILL_NO_SEMICOLON, n, *this->currentPos);
					}
					switch (n) {
					case '(':
						if (this->prevTokenState == NWK_LABEL) {
							throw OTCParsingError(_ILL_AFTER_LABEL, n, *this->currentPos);
						}
						if (this->prevTokenState == NWK_BRANCH_INFO) {
							throw OTCParsingError(_ILL_AFTER_BRANCH_INFO, n, *this->currentPos);
						}
						if (this->prevTokenState == NWK_CLOSE) {
							throw OTCParsingError(_ILL_AFTER_CLOSE, n, *this->currentPos);
						}
						if (this->prevTokenState == NWK_COLON) {
							throw OTCParsingError(_ILL_AFTER_COLON, n, *this->currentPos);
						}
						this->numUnclosedParens += 1;
						this->currTokenState = NWK_OPEN;
						this->currWord.assign(1, '(');
						return;
					case ')':
						if (this->prevTokenState == NWK_OPEN) {
							throw OTCParsingError(_ILL_AFTER_OPEN, n, *this->currentPos);
						}
						if (this->prevTokenState == NWK_COMMA) {
							throw OTCParsingError(_ILL_AFTER_COMMA, n, *this->currentPos);
						}
						if (this->prevTokenState == NWK_COLON) {
							throw OTCParsingError(_ILL_AFTER_COLON, n, *this->currentPos);
						}
						if (this->numUnclosedParens <= 0) {
							assert(this->numUnclosedParens == 0);
							throw OTCParsingError("Too many close parentheses", n, *this->currentPos);
						}
						this->numUnclosedParens -= 1;
						this->currTokenState = NWK_CLOSE;
						this->currWord.assign(1, ')');
						return;
					case ',':
						if (this->prevTokenState == NWK_OPEN) {
							throw OTCParsingError(_ILL_AFTER_OPEN, n, *this->currentPos);
						}
						if (this->prevTokenState == NWK_COLON) {
							throw OTCParsingError(_ILL_AFTER_COLON, n, *this->currentPos);
						}
						this->currTokenState = NWK_COMMA;
						this->currWord.assign(1, ':');
						return;
					case ':':
						if (this->prevTokenState == NWK_LABEL || this->prevTokenState == NWK_CLOSE) {
							this->currTokenState = NWK_COMMA;
							this->currWord.assign(1, ':');
							return;
						}
						if (this->prevTokenState == NWK_BRANCH_INFO) {
							throw OTCParsingError(_ILL_AFTER_BRANCH_INFO, n, *this->currentPos);
						}
						this->throwSCCErr(n);
					case ';':
						if (this->numUnclosedParens != 0) {
							throw OTCParsingError(_EARLY_SEMICOLON, n, *this->currentPos);
						}
						if (this->prevTokenState == NWK_COLON) {
							throw OTCParsingError(_ILL_AFTER_COLON, n, *this->currentPos);
						}
						if (this->prevTokenState == NWK_LABEL || this->prevTokenState == NWK_CLOSE) {
							this->currTokenState = NWK_NOT_IN_TREE;
							this->currWord.assign(1, ';');
							return;
						}
						this->throwSCCErr(n);
					case '[':
						this->finishReadingComment();
						break;
					case '\'':
						this->finishReadingQuotedStr();
						return;
					default:
						assert(false);
						throw OTCError("should be unreachable");
					}
				}
			}
		}
	}
}

} //namespace otc
