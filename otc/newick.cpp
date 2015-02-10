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
					this->finishReadingUnquoted();
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
							throw OTCParsingError(_ILL_AFTER_OPEN, n, *this->currentPos);
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
						break;
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
						break;
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