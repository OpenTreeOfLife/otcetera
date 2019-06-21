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
    uint64_t top, mid, bot;
    static constexpr unsigned int END_LETTER_INDEX = LETTER_INDEX_OF_FIRST_BIT_IN_THIRD_WORD + 64 - NUM_INDEX_BITS;

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
    uint64_t top, bot;
    static constexpr unsigned int END_LETTER_INDEX = LETTER_INDEX_OF_FIRST_BIT_IN_THIRD_WORD - NUM_INDEX_BITS;
    
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
    uint64_t top;
    static constexpr unsigned int END_LETTER_INDEX = LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD - NUM_INDEX_BITS;
    
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

//@TODO use lookup table
inline void fill_letter_and_node_indices(unsigned char curr_byte,
                                         int offset,
                                         vec_ind_pair_t & ret,
                                         uint64_t & node_index) {
    if (curr_byte == 0) {
        return;
    }                             
    unsigned char curr_bit = FIRST_BIT_OF_BYTE;
    for (unsigned char i = 0; i < 8; ++i) {
        //std::cerr << "fill_letter_and_node_indices byte=" << std::hex << (unsigned int)curr_byte << " bit=" << std::hex << (unsigned int)curr_bit << '\n' << std::dec;
        if (curr_byte & curr_bit) {
            ret.push_back(ind_pair_t{i + offset, node_index++});
        }
        curr_bit >>= 1;
    }
}

inline void fill_letter_and_node_indices_64(uint64_t masked_workspace,
                                            int offset,
                                            vec_ind_pair_t & ret,
                                            uint64_t & node_index) {
    int bitshift = 56;
    uint64_t blot_out_masker;
    for (unsigned char i = 0U; i < 8; ++i) {
        if (masked_workspace == 0) {
            return;
        }
        unsigned char curr_byte = masked_workspace >> bitshift;
        if (curr_byte != 0) {
            blot_out_masker = curr_byte;
            blot_out_masker <<= bitshift;
            // 0 out the bits in masked_workspace that we're dealing with in this iteration.
            masked_workspace ^= blot_out_masker; 
            fill_letter_and_node_indices(curr_byte, offset, ret, node_index);
        }
        bitshift -= 8;
        offset += 8;
    }
}

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
inline void CTrieNode<CTrie3NodeData>::flag_letter(unsigned int i) {
    uint64_t bit = ONE_64;
    //log_state();
    if (i < LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD) {
        const uint64_t shifted = (bit << (LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD - 1 - i));
        data.top |= shifted;
    } else if (i < LETTER_INDEX_OF_FIRST_BIT_IN_THIRD_WORD) {
        bit <<= (LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD - 1 - i);
        data.mid |= bit;
    } else {
        assert(i < DATA_TYPE::END_LETTER_INDEX);
        bit <<= (LETTER_INDEX_OF_FIRST_BIT_IN_THIRD_WORD -1 - i);
        data.bot |= bit;
    }
}

template <>
inline vec_ind_pair_t CTrieNode<CTrie3NodeData>::get_letter_and_node_indices_for_on_bits() const {
    assert(!is_terminal());
    vec_ind_pair_t ret;
    ret.reserve(DATA_TYPE::END_LETTER_INDEX);
    uint64_t node_index = get_index();
    uint64_t masked = data.top & TOP_LETTER_MASK;
    fill_letter_and_node_indices_64(masked, LETTER_INDEX_OF_FIRST_BIT_IN_FIRST_WORD, ret, node_index);
    masked = data.mid;
    fill_letter_and_node_indices_64(masked, LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD, ret, node_index);
    masked = data.bot & BOTTOM_LETTER_MASK;
    fill_letter_and_node_indices_64(masked, LETTER_INDEX_OF_FIRST_BIT_IN_THIRD_WORD, ret, node_index);
    return ret;
}

template <>
inline void CTrieNode<CTrie2NodeData>::flag_letter(unsigned int i) {
    uint64_t bit = ONE_64;
    if (i < LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD) {
        bit <<= (LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD - 1 - i);
        data.top |= bit;
    }  else {
        assert(i < DATA_TYPE::END_LETTER_INDEX);
        bit <<= (LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD - 1 - i);
        data.bot |= bit;
    }
}

template <>
inline vec_ind_pair_t CTrieNode<CTrie2NodeData>::get_letter_and_node_indices_for_on_bits() const {
    //std::cerr << "get_letter_and_node_indices_for_on_bits top="
    //          << std::hex << top << " bot=" << std::hex << bot << std::dec << '\n';
    assert(!is_terminal());
    vec_ind_pair_t ret;
    ret.reserve(DATA_TYPE::END_LETTER_INDEX);
    u_int64_t node_index = get_index();
    uint64_t masked = data.top & TOP_LETTER_MASK;
    fill_letter_and_node_indices_64(masked, LETTER_INDEX_OF_FIRST_BIT_IN_FIRST_WORD, ret, node_index);
    masked = data.bot & BOTTOM_LETTER_MASK;
    fill_letter_and_node_indices_64(masked, LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD, ret, node_index);
    return ret;
}


template <>
inline void CTrieNode<CTrie1NodeData>::flag_letter(unsigned int i) {
    uint64_t bit = ONE_64;
    assert(i < DATA_TYPE::END_LETTER_INDEX);
    bit <<= (LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD - 1 - i);
    data.top |= bit;
} 

template <>
inline vec_ind_pair_t CTrieNode<CTrie1NodeData>::get_letter_and_node_indices_for_on_bits() const {
    //std::cerr << "get_letter_and_node_indices_for_on_bits top="
    //          << std::hex << top << " bot=" << std::hex << bot << std::dec << '\n';
    assert(!is_terminal());
    vec_ind_pair_t ret;
    ret.reserve(DATA_TYPE::END_LETTER_INDEX);
    u_int64_t node_index = get_index();
    uint64_t masked = data.top & TOP_LETTER_MASK;
    fill_letter_and_node_indices_64(masked, LETTER_INDEX_OF_FIRST_BIT_IN_FIRST_WORD, ret, node_index);
    return ret;
}

using CTrie3Node = CTrieNode<CTrie3NodeData>;
using CTrie2Node = CTrieNode<CTrie2NodeData>;
using CTrie1Node = CTrieNode<CTrie1NodeData>;

} // namespace otc
#endif
