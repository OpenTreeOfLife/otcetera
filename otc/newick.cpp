#include "otc/newick.h"
#include <cstring>
#include <string>
namespace otc {

static const char * _EARLY_SEMICOLON = "Unexpected ; with open parentheses not balanced.";
static const char * _ILL_AFTER_CLOSE = "Illegal character after \")\" character. Expecting \",\" or a label or a colon.";
static const char * _ILL_AFTER_LABEL = "Illegal character after label. Expecting ( or a label.";
static const char * _ILL_AFTER_BRANCH_INFO = "Illegal character after branch info. Expecting ( or a label.";
static const char * _ILL_AFTER_OPEN = "Illegal character after \"(\" character. Expecting ( or a label.";
static const char * _ILL_AFTER_COMMA = "Illegal character after \",\" character. Expecting ( or a label.";
static const char * _ILL_AFTER_COLON = "Illegal character after \":\" character. Expecting a branch length.";
static const char * _ILL_NO_SEMICOLON = "Expecting ; after a newick description.";
static const char * _ILL_FIRST_CHAR = "Expecting a newick tree to start with \"(\".";

void NewickTokenizer::iterator::on_label_exit(char n, bool fromWS) {
    bool whitespaceFound = fromWS;
    if (std::strchr("(),:;", n) == nullptr) {
        if (!std::isgraph(n)) {
            whitespaceFound = true;
            if (!advance_to_next_non_whitespace(n)) {
                return;
            }
            if (std::strchr("(),:;", n) != nullptr) { // exit stripping whitespace
                this->push(n);
                return;
            }
            this->current_word += ' '; // replace any series of one white space with one " ". Not sure what to do, this seems reasonable.
        }
        if (n == '[') {
            this->advance_reader_one_logical_char(n);
            assert(n == '[');
            finish_reading_comment();
            n = this->peek();
            on_label_exit(n, false);
            return;
        }
        LOG(WARNING) << "Unexpected continuation of a label after a quoted string in newick. Next character is: \"" << n << "\"";
        if (whitespaceFound) {
            this->current_word += ' ';
        }
        if (n == '\'') {
            this->advance_reader_one_logical_char(n);
            assert(n == '\'');
            finish_reading_quoted_str();
            return;
        }
        this->current_word += n;
        this->advance_reader_one_logical_char(n);
        this->finish_reading_unquoted(true);
        return;
    }
}
void NewickTokenizer::iterator::finish_reading_quoted_str(){
    for (;;) {
        char c;
        if(!advance_reader_one_logical_char(c)) {
            throw OTCParsingError("Unexpected EOF in quoted string", '\0', (*this->current_pos));
        }
        if (c == '\'') {
            char n;
            if (!advance_reader_one_logical_char(n)) {
                    throw OTCParsingError("Unexpected EOF at the end of a quoted string. Expecting ; if this is the end of the tree.", '\0', (*this->current_pos));
            }
            if (n == '\'') {
                this->current_word += c;
            } else {
                this->push(n);
                this->on_label_exit(n, false);
                return;
            }
        } else {
            this->current_word += c;
        }
    }
}
void NewickTokenizer::iterator::finish_reading_unquoted(bool continuingLabel){
    for (;;) {
        char c;
        if(!advance_reader_one_logical_char(c)) {
            throw OTCParsingError("Unexpected EOF in label. Expecting a ; to end a newick.", '\0', (*this->current_pos));
        }
        if (!isgraph(c)) {
            if (!advance_to_next_non_whitespace(c)) {
                throw OTCParsingError("Unexpected EOF in label. Expecting a ; to end a newick.", '\0', (*this->current_pos));
            }
            if (continuingLabel) {
                LOG(WARNING) << "Whitespace found in unquoted label - will be converted to a single space.";
            }
            this->push(c);
            this->on_label_exit(c, true);
            return;
        }
        if (std::strchr("(),:;", c) != nullptr) {
            this->push(c);
            return;
        }
        if (c == '\'') {
            LOG(WARNING) << "single-quoted string found in unquoted label.";
            this->finish_reading_quoted_str();
            return;
        } else if (c == '[') {
            finish_reading_comment();
        } else if (c == '_') {
            this->current_word += ' ';
        } else {
            this->current_word += c;
        }
    }
}
void NewickTokenizer::iterator::finish_reading_comment(){
    auto numOpenComments = 1U;
    std::string comment;
    for (;;) {
        char c;
        if(!advance_reader_one_logical_char(c)) {
            throw OTCParsingError("Unexpected EOF in comment", '\0', (*this->current_pos));
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
bool NewickTokenizer::iterator::advance_to_next_non_whitespace(char & c){
    for (;;) {
        if(!advance_reader_one_logical_char(c)) {
            return false;
        }
        if (std::isgraph(c)) {
            return true;
        }
    }
}

void NewickTokenizer::iterator::throw_scc_err(char n) const {
    if (this->previous_token_state == NWK_OPEN) {
        throw OTCParsingError(_ILL_AFTER_OPEN, n, *this->current_pos);
    }
    if (this->previous_token_state == NWK_COLON) {
        throw OTCParsingError(_ILL_AFTER_COLON, n, *this->current_pos);
    }
    if (this->previous_token_state == NWK_COMMA) {
        throw OTCParsingError(_ILL_AFTER_OPEN, n, *this->current_pos);
    }
    assert(false);
    throw OTCError("should be unreachable");
}

void NewickTokenizer::iterator::consume_next_token() {
    char n;
    // a loop to accumulate comments
    while (true) {
        if (!advance_to_next_non_whitespace(n)) {
            LOG(TRACE) << "in false consume_next_token with n =\"" << n << "\"";
            if (this->current_token_state == NWK_NOT_IN_TREE || this->current_token_state == NWK_SEMICOLON) {
                return;
            }
            if (this->num_unclosed_parens > 0) {
                throw OTCParsingError("Unexpected EOF while tree was still being read more open parentheses than closed parentheses read.", '\0', *this->current_pos);
            }
            throw OTCParsingError("Unexpected EOF. Semicolon expected at the end of the newick.", '\0', *this->current_pos);
        } else {
            LOG(TRACE) << "in consume_next_token with n =\"" << n << "\"";
            if (this->previous_token_state == NWK_NOT_IN_TREE) {
                if (n == '(') {
                    this->current_token_state = NWK_OPEN;
                    this->num_unclosed_parens += 1;
                    this->current_word.assign(1, '(');
                    return;
                } else if (n == '[') {
                    this->finish_reading_comment();
                } else {
                    throw OTCParsingError(_ILL_FIRST_CHAR, n, *this->current_pos);
                }
            }
            else {
                if (std::strchr("(),:;[\'", n) == nullptr) {
                    this->current_token_state = (this->previous_token_state == NWK_COLON ? NWK_BRANCH_INFO : NWK_LABEL);
                    this->current_word.assign(1, n);
                    this->finish_reading_unquoted(false);
                    return;
                } else {
                    switch (n) {
                    case '(':
                        if (this->previous_token_state == NWK_LABEL) {
                            throw OTCParsingError(_ILL_AFTER_LABEL, n, *this->current_pos);
                        }
                        if (this->previous_token_state == NWK_BRANCH_INFO) {
                            throw OTCParsingError(_ILL_AFTER_BRANCH_INFO, n, *this->current_pos);
                        }
                        if (this->previous_token_state == NWK_CLOSE) {
                            throw OTCParsingError(_ILL_AFTER_CLOSE, n, *this->current_pos);
                        }
                        if (this->previous_token_state == NWK_COLON) {
                            throw OTCParsingError(_ILL_AFTER_COLON, n, *this->current_pos);
                        }
                        if (this->num_unclosed_parens <= 0) {
                            assert(this->num_unclosed_parens == 0);
                            throw OTCParsingError(_ILL_NO_SEMICOLON, n, *this->current_pos);
                        }
                        this->num_unclosed_parens += 1;
                        this->current_token_state = NWK_OPEN;
                        this->current_word.assign(1, '(');
                        return;
                    case ')':
                        if (this->previous_token_state == NWK_OPEN) {
                            throw OTCParsingError(_ILL_AFTER_OPEN, n, *this->current_pos);
                        }
                        if (this->previous_token_state == NWK_COMMA) {
                            throw OTCParsingError(_ILL_AFTER_COMMA, n, *this->current_pos);
                        }
                        if (this->previous_token_state == NWK_COLON) {
                            throw OTCParsingError(_ILL_AFTER_COLON, n, *this->current_pos);
                        }
                        if (this->num_unclosed_parens <= 0) {
                            assert(this->num_unclosed_parens == 0);
                            throw OTCParsingError("Too many close parentheses", n, *this->current_pos);
                        }
                        this->num_unclosed_parens -= 1;
                        this->current_token_state = NWK_CLOSE;
                        this->current_word.assign(1, ')');
                        return;
                    case ',':
                        if (this->previous_token_state == NWK_OPEN) {
                            throw OTCParsingError(_ILL_AFTER_OPEN, n, *this->current_pos);
                        }
                        if (this->previous_token_state == NWK_COLON) {
                            throw OTCParsingError(_ILL_AFTER_COLON, n, *this->current_pos);
                        }
                        if (this->num_unclosed_parens <= 0) {
                            assert(this->num_unclosed_parens == 0);
                            throw OTCParsingError(_ILL_NO_SEMICOLON, n, *this->current_pos);
                        }
                        if (this->previous_token_state == NWK_COMMA) {
                            throw OTCParsingError(_ILL_AFTER_COMMA, n, *this->current_pos);
                        }
                        this->current_token_state = NWK_COMMA;
                        this->current_word.assign(1, ',');
                        return;
                    case ':':
                        if (this->previous_token_state == NWK_LABEL || this->previous_token_state == NWK_CLOSE) {
                            this->current_token_state = NWK_COLON;
                            this->current_word.assign(1, ':');
                            return;
                        }
                        if (this->previous_token_state == NWK_BRANCH_INFO) {
                            throw OTCParsingError(_ILL_AFTER_BRANCH_INFO, n, *this->current_pos);
                        }
                        this->throw_scc_err(n);
                    case ';':
                        if (this->num_unclosed_parens != 0) {
                            throw OTCParsingError(_EARLY_SEMICOLON, n, *this->current_pos);
                        }
                        if (this->previous_token_state == NWK_COLON) {
                            throw OTCParsingError(_ILL_AFTER_COLON, n, *this->current_pos);
                        }
                        if (this->previous_token_state == NWK_LABEL || this->previous_token_state == NWK_BRANCH_INFO || this->previous_token_state == NWK_CLOSE) {
                            this->current_token_state = NWK_SEMICOLON;
                            this->current_word.assign(1, ';');
                            return;
                        }
                        this->throw_scc_err(n);
                    case '[':
                        this->finish_reading_comment();
                        break;
                    case '\'':
                        // setting of the current_token_state here assumes that the quoted token will
                        //      be a label not a part of the newick tree syntax.
                        // One could argue that we should allow all tokens to be quoted. 
                        //      in practice, no one seems to quote their () , : or ; characters.
                        this->current_token_state = (this->previous_token_state == NWK_COLON ? NWK_BRANCH_INFO : NWK_LABEL);
                        this->finish_reading_quoted_str();
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
