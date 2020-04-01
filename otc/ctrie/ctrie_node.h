#ifndef OTC_CTRIE_NODE_H
#define OTC_CTRIE_NODE_H
#include <cstddef>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <bitset>
#include "otc/error.h"
#include "otc/otc_base_includes.h"

namespace otc {

using stored_index_t = unsigned char;

// The CTrieXNode classes are the elements stored in a vector by the CTrie class.
// To allow for random access into the vector, each node is the same size for 
//  a type, despite the fact that there are 2 distinct types of nodes:
//      1. Terminal nodes just hold an index to the suffix for the trie.
//      2. Internal nodes hold flags for what letters are next in the trie and
//           the index of the first daughter.
// The first bit in the data field indicates whether or not the node is a terminal.
//
// highest bit: is terminal node
//  If the top bit is is 0:
//      second highest bit: has key that terminates with this node
//      letter bits = the X bits until the bottom uint64_t.
//  If the top bit is is 1:
//      second highest and letter bits are unused.
//      If top bit is 1, then 
// The lowest NUM_INDEX_BITS of the bottom uint64_t are used
//   to hold an index. Either:
//      1. The first of the next daughter nodes 
//         (in the ctrie vector of nodes) if top bit=0, or
//      2. the index of the offset into the ctrie suffix 
//          vector that holds the suffix for the string. 
constexpr unsigned char NUM_INDEX_BITS = 32;
constexpr uint64_t ZERO_64 = 0;
constexpr uint64_t ONE_64 = 1;
constexpr uint64_t HIGHEST_BIT = ONE_64 << 63;
constexpr uint64_t SECOND_HIGHEST_BIT = ONE_64 << 62;
constexpr uint64_t TOP_LETTER_MASK = SECOND_HIGHEST_BIT - 1;
constexpr uint64_t INDEX_MASK = (ONE_64 << NUM_INDEX_BITS) - ONE_64;
constexpr uint64_t COMP_INDEX_MASK = ~INDEX_MASK;
constexpr uint64_t BOTTOM_LETTER_MASK = COMP_INDEX_MASK;
constexpr stored_index_t NO_MATCHING_CHAR_CODE = 255;
constexpr unsigned char top_first_byte = 0x03F; // 64 and 128 bit used for top 2 flags
constexpr unsigned char FULL_BYTE = 0x0FF;
constexpr unsigned char FIRST_BIT_OF_BYTE = 0x080;
constexpr int LETTER_INDEX_OF_FIRST_BIT_IN_FIRST_WORD = -2; // 2 bits for flags
constexpr unsigned int LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD = 64 + LETTER_INDEX_OF_FIRST_BIT_IN_FIRST_WORD;
constexpr unsigned int LETTER_INDEX_OF_FIRST_BIT_IN_THIRD_WORD = 64 + LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD;

class CTrie3NodeData {
public:
    static constexpr unsigned int END_LETTER_INDEX = LETTER_INDEX_OF_FIRST_BIT_IN_THIRD_WORD + 64 - NUM_INDEX_BITS;

    uint64_t top, mid, bot;

    CTrie3NodeData() :top{ZERO_64}, mid{ZERO_64}, bot{ZERO_64} {
    }
    uint64_t & get_flag_word() {
        return top;
    }
    const uint64_t & get_flag_word_const() const {
        return top;
    }
    uint64_t & get_middle_word() {
        return mid;
    }
    const uint64_t & get_middle_word_const() const {
        return mid;
    }
    uint64_t & get_index_word() {
        return bot;
    }
    const uint64_t & get_index_word_const() const {
        return bot;
    }

    void db_write_state(std::ostream &out) const {
       out << "top=" << std::bitset<64>{top} << " mid=" << std::bitset<64>{mid}  << " bot=" << std::bitset<64>{bot} ;
    }
};

class CTrie2NodeData {
public:
    static constexpr unsigned int END_LETTER_INDEX = LETTER_INDEX_OF_FIRST_BIT_IN_THIRD_WORD - NUM_INDEX_BITS;

    uint64_t top, bot;
    
    CTrie2NodeData() :top{ZERO_64}, bot{ZERO_64} {
    }
    uint64_t & get_flag_word() {
        return top;
    }
    const uint64_t & get_flag_word_const() const {
        return top;
    }
    uint64_t & get_index_word() {
        return bot;
    }
    const uint64_t & get_index_word_const() const {
        return bot;
    }

    void db_write_state(std::ostream &out) const {
       out << "top=" << std::bitset<64>{top} << " bot=" << std::bitset<64>{bot} ;
    }
};

class CTrie1NodeData {
public:
    static constexpr unsigned int END_LETTER_INDEX = LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD - NUM_INDEX_BITS;

    uint64_t top;
    
    CTrie1NodeData() :top{ZERO_64} {
    }
    uint64_t & get_flag_word() {
        return top;
    }
    const uint64_t & get_flag_word_const() const {
        return top;
    }
    uint64_t & get_index_word() {
        return top;
    }
    const uint64_t & get_index_word_const() const {
        return top;
    }

    void db_write_state(std::ostream &out) const {
       out << "top=" << std::bitset<64>{top}  ;
    }
};

using ind_pair_t = std::pair<stored_index_t, uint64_t>;
using vec_ind_pair_t = std::vector<ind_pair_t>;

template <typename T>
class CTrieNode {
private:
    T data;
public:
    using DATA_TYPE = T;

    CTrieNode() {
    }

    uint64_t get_index() const {
        return INDEX_MASK & data.get_index_word_const();
    }
    
    void log_state() const {
       std::cerr << " CTrieNode( "; data.db_write_state(std::cerr); std::cerr << ")\n";
    }

    void flag_as_key_terminating() {
        data.get_flag_word() |= SECOND_HIGHEST_BIT;
    }

    bool is_key_terminating() const {
        return data.get_flag_word_const() & SECOND_HIGHEST_BIT;
    }

    void flag_as_terminal() {
        data.get_flag_word()  |= HIGHEST_BIT;
    }

    bool is_terminal() const {
        return data.get_flag_word_const() & HIGHEST_BIT;
    }
    
    void set_index(std::size_t index) {
        uint64_t ind = index;
        ind &= INDEX_MASK;
        if (ind != (uint64_t)index) {
            throw OTCError() << "not enough index field to hold pos = " << index;
        }
        auto & word = data.get_index_word();
        word &= COMP_INDEX_MASK; // sets to 0 any bits for the index
        word |= ind;
    }

    void flag_as_suffix(std::size_t pos) {
        flag_as_terminal();
        set_index(pos);
    }

    void set_first_child_index(std::size_t index) {
        assert(!is_terminal());
        set_index(index);
    }


    void flag_letter(unsigned int i);
    vec_ind_pair_t get_letter_and_node_indices_for_on_bits() const;
};


template <>
void CTrieNode<CTrie3NodeData>::flag_letter(unsigned int i);

template <>
vec_ind_pair_t CTrieNode<CTrie3NodeData>::get_letter_and_node_indices_for_on_bits() const;

template <>
void CTrieNode<CTrie2NodeData>::flag_letter(unsigned int i);

template <>
vec_ind_pair_t CTrieNode<CTrie2NodeData>::get_letter_and_node_indices_for_on_bits() const;

template <>
void CTrieNode<CTrie1NodeData>::flag_letter(unsigned int i);

template <>
vec_ind_pair_t CTrieNode<CTrie1NodeData>::get_letter_and_node_indices_for_on_bits() const;

using CTrie3Node = CTrieNode<CTrie3NodeData>;
using CTrie2Node = CTrieNode<CTrie2NodeData>;
using CTrie1Node = CTrieNode<CTrie1NodeData>;

} // namespace otc
#endif
