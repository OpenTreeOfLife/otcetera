#ifndef OTCETERA_NEWICK_TOKENIZER_H
#define OTCETERA_NEWICK_TOKENIZER_H
#include <iostream>
#include <fstream>
#include <stack>
#include <map>
#include <stdexcept>
#include "otc/otc_base_includes.h"
#include "otc/error.h"
#include "otc/util.h"

namespace otc {

struct ParsingRules {
    const OttIdSet * ott_id_validator = nullptr;
    bool include_internal_nodes_in_des_id_sets = false;
    bool set_ott_idForInternals = true; // TEMP. Do we want to just remove them after processing, rather than adding a conditional to the parsing logic...
    const std::map<OttId, OttId> * id_remapping = nullptr;
    bool prune_unrecognized_input_tips = false;
    bool require_ott_ids = true;  // Every label must include an OttId
    bool set_ott_ids = true;      // Read and set OttIds for labels that have them.
};

typedef std::shared_ptr<const std::string> ConstStrPtr;
struct FilePosStruct {
    FilePosStruct() = default;
    FilePosStruct(const FilePosStruct&) = default;
    FilePosStruct(const ConstStrPtr fp):filepath(fp) {}
    FilePosStruct(std::size_t position, std::size_t line, std::size_t column, const ConstStrPtr& fp)
        :pos(position),
        lineNumber(line),
        colNumber(column),
        filepath(fp) {
    }

    void set_location_in_file(const FilePosStruct & other) {
        pos = other.pos;
        lineNumber = other.lineNumber;
        colNumber = other.colNumber;
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

    std::size_t pos = 0;
    std::size_t lineNumber = 0;
    std::size_t colNumber = 0;
    ConstStrPtr filepath = nullptr;
};

class OTCParsingError: public OTCError {
    public:
        OTCParsingError(const char * msg)
	{
	    (*this)<<"parsing newick: "<<msg;
	}

        OTCParsingError(const char * msg, const FilePosStruct & position)
	    :pos(position)
	{
	    (*this)<<"parsing newick: "<<msg<<" "<<pos.describe();
	}

        OTCParsingError(const char * msg, char offending, const FilePosStruct & position)
            :OTCError(msg),
            frag(msg),
            pos(position) {
                if (offending != '\0') {
                    offendingStr.assign(1, offending);
                }
                message = this->generate_message();
            }
        OTCParsingError(const char * msg, const std::string & offending, const FilePosStruct & position)
            :OTCError(msg),
            frag(msg),
            offendingStr(offending),
            pos(position) {
                message = this->generate_message();
            }
        std::string generate_message() const noexcept {
            try {
                std::string m = "Error found \"";
                m += this->offendingStr;
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
        std::string offendingStr;
        const FilePosStruct pos;
};

class OTCParsingContentError: public OTCError {
    public:
        OTCParsingContentError(const char * msg, const FilePosStruct & position)
            :OTCError(msg),
            pos(position) {
                message = this->generate_message(msg);
            }
        std::string generate_message(const char *frag) const noexcept {
            try {
                std::string m = "Error found: ";
                m += frag;
                m += " ";
                m += this->pos.describe();
                return m;
            } catch (...) {// to guarantee noexcept...
            }
            try {
                return std::string{frag};
            } catch (...) {
                return std::string{};
            }
        }
    private:
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
        NewickTokenizer(std::istream &inp, const FilePosStruct & initialPos)
            :input_stream(inp),
            initPos(initialPos) {
        }
        class iterator;
        class Token {
            public:
                const std::string & content() const {
                    return this->token_content;
                }
                const std::vector<std::string> & comment_vec() const {
                    return this->comments;
                }
                const FilePosStruct & get_start_pos() const {
                    return this->start_pos;
                }
            private:
                Token(const std::string &content,
                      const FilePosStruct & startPosition,
                      const FilePosStruct & endPosition,
                      const std::vector<std::string> &embeddedComments,
                      newick_token_state_t tokenState)
                    :token_content(content),
                    start_pos(startPosition),
                    end_pos(endPosition),
                    comments(embeddedComments),
                    state(tokenState) {
                    //LOG(TRACE) << "created token for \"" << content << "\"";
                }
            public:
                const std::string token_content;
                const FilePosStruct start_pos;
                const FilePosStruct end_pos;
                const std::vector<std::string> comments;
                const newick_token_state_t state;

                friend class NewickTokenizer::iterator;
        };

        class iterator : std::forward_iterator_tag {
            public:
                bool operator==(const iterator & other) const {
                    //LOG(TRACE) << "Equality test at_end = " << this->at_end << " other.at_end = " << other.at_end;
                    if (other.at_end) {
                        return this->at_end;
                    }
                    if (this->at_end) {
                        return false;
                    }
                    return (this->current_pos == other.current_pos) && (&(this->input_stream) == &(other.input_stream));
                }
                bool operator!=(const iterator & other) const {
                    //LOG(TRACE) << "Inequality test at_end = " << this->at_end << " other.at_end = " << other.at_end << '\n';
                    return !(*this == other);
                }
                Token operator*() const {
                    //LOG(TRACE) << "* operator";
                    return Token(current_word, *prev_pos, *current_pos, comments, current_token_state);
                }
                iterator & operator++() {
                    //LOG(TRACE) << "increment";
                    if (this->at_end) {
                        throw std::out_of_range("Incremented a dead NewickTokenizer::iterator");
                    }
                    if (this->current_token_state == NWK_SEMICOLON) {
                        this->at_end = true;
                        this->reset_token();
                    } else {
                        this->reset_token();
                        this->consume_next_token();
                    }
                    return *this;
                }
                const FilePosStruct & get_curr_pos() const {
                    return *this->current_pos;
                }
            private:
                void consume_next_token();
                bool advance_to_next_non_whitespace(char &);
                void finish_reading_comment();
                void finish_reading_unquoted(bool continuingLabel);
                void finish_reading_quoted_str();
                void on_label_exit(char nextChar, bool enteringFromWhitespace);
                char peek() {
                    if (!pushed.empty()) {
                        return pushed.top();
                    }
                    char c = static_cast<char>(this->input_stream.rdbuf()->sgetc());
                    pushed.push(c);
                    return c;
                }
                void push(char c) {
                    this->pushed.push(c);
                    if (c == '\n') {
                        this->current_pos->colNumber = last_line_ind;
                        this->current_pos->lineNumber -= 1;
                    } else {
                        this->current_pos->colNumber -= 1;
                    }
                    this->current_pos->pos -= 1;
                }
                void throw_scc_err(char c) const __attribute__ ((noreturn));
                //deals with \r\n as \n Hence "LogicalChar"
                bool advance_reader_one_logical_char(char & c) {
                    if (!pushed.empty()) {
                        c = pushed.top();
                        pushed.pop();
                    } else {
                        if (this->at_end) {
                            c = EOF;
                            return false;
                        }
                        c = static_cast<char>((this->input_stream.rdbuf())->sbumpc());
                    }
                    if (c == EOF) {
                        this->at_end = true;
                    } else {
                        this->current_pos->pos += 1;
                        if (13 == c || 10 == c) {
                            if (13 == c ) { // deal with \r\n as a newline
                                if (this->input_stream.rdbuf()->sgetc() == 10) {//peeks at the next char
                                    (input_stream.rdbuf())->sbumpc();
                                    this->current_pos->pos += 1;
                                }
                            }
                            c = '\n';
                            this->last_line_ind = this->current_pos->colNumber;
                            this->current_pos->colNumber = 0;
                            this->current_pos->lineNumber += 1;
                        } else {
                            this->current_pos->colNumber += 1;
                        }
                    }
                    return !this->at_end;
                }
                void reset_token() {
                    this->current_word.clear();
                    this->previous_token_state = this->current_token_state;
                    std::swap(current_pos, prev_pos);
                    current_pos->pos = prev_pos->pos;
                    current_pos->lineNumber = prev_pos->lineNumber;
                    current_pos->colNumber = prev_pos->colNumber;
                    comments.clear();
                }
                iterator(std::istream &inp, const FilePosStruct & initialPos)
                    :input_stream(inp),
                    input_filepath(initialPos.filepath),
                    at_end(!inp.good()),
                    first_pos_slot(initialPos),
                    second_pos_slot(initialPos.filepath),
                    current_token_state(NWK_NOT_IN_TREE),
                    previous_token_state(NWK_NOT_IN_TREE),
                    num_unclosed_parens(0) {
                    current_pos = &first_pos_slot;
                    prev_pos = &second_pos_slot;
                    //LOG(TRACE) << "create live";
                    ++(*this);
                }
                iterator(std::istream &inp) // USE in end() ONLY!
                    :input_stream(inp),
                    input_filepath(nullptr),
                    at_end(true),
                    first_pos_slot(nullptr),
                    second_pos_slot(nullptr),
                    current_token_state(NWK_NOT_IN_TREE),
                    previous_token_state(NWK_NOT_IN_TREE),
                    num_unclosed_parens(0) {
                    //LOG(TRACE) << "create dead";

                }
                std::istream & input_stream;
                ConstStrPtr input_filepath;
                bool at_end;
                FilePosStruct first_pos_slot;
                FilePosStruct second_pos_slot;
                FilePosStruct * current_pos; // alias
                FilePosStruct * prev_pos;    // alias
                std::string current_word;
                newick_token_state_t current_token_state;
                newick_token_state_t previous_token_state;
                long num_unclosed_parens;
                std::stack<char> pushed;
                std::vector<std::string> comments;
                std::size_t last_line_ind;
                friend class NewickTokenizer;
        };
        iterator begin() {
            iterator b(this->input_stream, initPos);
            return b;
        }
        iterator end() {
            //LOG(TRACE) << "NewickTokenizer.end()";
            return iterator(this->input_stream);
        }
    private:
        std::istream & input_stream;
        FilePosStruct initPos;
};

/// Get OTT Id from a string of the form (ott######) or (.......[ \t_]ott#####).
inline long long_ott_id_from_name(const std::string & n) {
    if (n.empty()) {
        return -1;
    }
    const unsigned long lastInd = n.length() - 1;
    unsigned long currInd = lastInd;
    const char * c = n.c_str();
    if (strchr("0123456789", c[currInd]) == 0) {
        return -2;
    }
    while (currInd > 0) {
        --currInd;
        if (strchr("0123456789", c[currInd]) == 0) {
            ++currInd;
            break;
        }
    }
    if (currInd < 3) {
        return -2;
    }
    if (currInd >= 3 and strncmp(c+currInd-3,"ott",3) != 0) {
        return -2;
    }
    // Valid separators between ott####### and previous characters.
    if (currInd > 3 and strchr("_ \t",c[currInd-4]) == 0) {
        return -2;
    }
    long conv = -2;
    auto r = char_ptr_to_long(c + currInd, &conv);
    assert(r);
    return conv;
}

/// Get the OTT Id from a string ###### consisting of digits only, with no whitespace or other characters.
inline long string_to_long_ott_id(const std::string & n) {
    if (n.empty()) {
        return -1;
    }
    long conv = -2;
    // If conversion does not consume the entire string, conv is unchanged.
    char_ptr_to_long(n.c_str(), &conv);
    return conv;
}

}// namespace otc
#endif
