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
constexpr uint8_t ZER0_8 = 0;
constexpr uint32_t ZERO_32 = 0;
constexpr uint64_t ZERO_64 = 0;
constexpr uint8_t ONE_8 = 1;
constexpr uint8_t HIGHEST_BIT8 = ONE_8 << 7;
constexpr uint8_t SECOND_HIGHEST_BIT8 = ONE_8 << 6;
constexpr uint8_t TOP_LETTER_MASK8 = SECOND_HIGHEST_BIT8 - 1;

constexpr stored_index_t NO_MATCHING_CHAR_CODE = 255;
constexpr int LETTER_INDEX_OF_FIRST_BIT_IN_FIRST_WORD = -2; // 2 bits for flags
constexpr unsigned char FIRST_BIT_OF_BYTE = 0x080;
/*
constexpr uint64_t ONE_64 = 1;
constexpr uint64_t HIGHEST_BIT = ONE_64 << 63;
constexpr uint64_t SECOND_HIGHEST_BIT = ONE_64 << 62;
constexpr uint64_t TOP_LETTER_MASK = SECOND_HIGHEST_BIT - 1;
constexpr uint64_t INDEX_MASK = (ONE_64 << NUM_INDEX_BITS) - ONE_64;
constexpr uint64_t COMP_INDEX_MASK = ~INDEX_MASK;
constexpr uint64_t BOTTOM_LETTER_MASK = COMP_INDEX_MASK;
constexpr unsigned char top_first_byte = 0x03F; // 64 and 128 bit used for top 2 flags
constexpr unsigned char FULL_BYTE = 0x0FF;
constexpr unsigned int LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD = 64 + LETTER_INDEX_OF_FIRST_BIT_IN_FIRST_WORD;
constexpr unsigned int LETTER_INDEX_OF_FIRST_BIT_IN_THIRD_WORD = 64 + LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD;
*/
class CTrie80NodeData {
    public:
    uint8_t letters[6];
    uint32_t index;
    static constexpr unsigned char NUM_LETTER_BYTES = 6;
    
    CTrie80NodeData()
        :letters{ZER0_8, ZER0_8, ZER0_8, ZER0_8, ZER0_8, ZER0_8},
        index{ZERO_32} {
    }
    uint8_t & get_flag_word() {
        return letters[0];
    }
    const uint8_t & get_flag_word_const() const {
        return letters[0];
    }
    uint32_t & get_index_word() {
        return index;
    }
    const uint32_t & get_index_word_const() const {
        return index;
    }

    void db_write_state(std::ostream &out) const {
       out << "top=";
       for (auto i = 0U; i < NUM_LETTER_BYTES; ++i) {
           out << std::bitset<8>{letters[i]};
       }
       out << " index=" << std::bitset<32>{index} ;
    }
};


class CTrie128NodeData {
    public:
    uint8_t letters[12];
    uint32_t index;
    static constexpr unsigned char NUM_LETTER_BYTES = 12;
    
    CTrie128NodeData()
        :letters{ZER0_8, ZER0_8, ZER0_8, ZER0_8, ZER0_8, ZER0_8, ZER0_8, ZER0_8, ZER0_8, ZER0_8, ZER0_8, ZER0_8},
        index{ZERO_32} {
    }
    uint8_t & get_flag_word() {
        return letters[0];
    }
    const uint8_t & get_flag_word_const() const {
        return letters[0];
    }
    uint32_t & get_index_word() {
        return index;
    }
    const uint32_t & get_index_word_const() const {
        return index;
    }

    void db_write_state(std::ostream &out) const {
       out << "top=";
       for (auto i = 0U; i < NUM_LETTER_BYTES; ++i) {
           out << std::bitset<8>{letters[i]};
       }
       out << " index=" << std::bitset<32>{index} ;
    }
};

using ind_pair_t = std::pair<stored_index_t, uint32_t>;
using vec_ind_pair_t = std::vector<ind_pair_t>;

//@TODO use lookup table
inline void fill_letter_and_node_indices(unsigned char curr_byte,
                                         int offset,
                                         vec_ind_pair_t & ret,
                                         uint32_t & node_index) {
    if (curr_byte == 0) {
        return;
    }                             
    unsigned char curr_bit = FIRST_BIT_OF_BYTE;
    for (unsigned char i = 0; i < 8; ++i) {
        if (curr_byte & curr_bit) {
            ret.push_back(ind_pair_t{i + offset, node_index++});
        }
        curr_bit >>= 1;
    }
}

template <typename T>
class CTrieNode {
    private:
    T data;
    public:
    static constexpr unsigned int MAX_LETTER_CAPACITY = 8*T::NUM_LETTER_BYTES - 2;
    using DATA_TYPE = T;
    CTrieNode() {
    }

    std::uint32_t get_index() const {
        return data.get_index_word_const();
    }
    
    void log_state(std::ostream & out) const {
        out << " CTrieNode( "; data.db_write_state(std::cerr); std::cerr << ")\n";
    }

    void flag_as_key_terminating() {
        data.get_flag_word() |= SECOND_HIGHEST_BIT8;
    }

    bool is_key_terminating() const {
        return data.get_flag_word_const() & SECOND_HIGHEST_BIT8;
    }

    void flag_as_terminal() {
        data.get_flag_word()  |= HIGHEST_BIT8;
    }

    bool is_terminal() const {
        return data.get_flag_word_const() & HIGHEST_BIT8;
    }
    
    void set_index(std::size_t index) {
        uint64_t ind = index;
        if (ind > UINT32_MAX) {
            throw OTCError() << "not enough index field to hold pos = " << index;
        }
        auto & word = data.get_index_word();
        word = (std::uint32_t) index;
    }

    void flag_as_suffix(std::size_t pos) {
        flag_as_terminal();
        set_index(pos);
    }

    void set_first_child_index(std::size_t index) {
        assert(!is_terminal());
        set_index(index);
    }

    void flag_letter(unsigned int i)  {
        unsigned int j = i + 2;
        auto byte_index = j / 8;
        auto bit_left_index = j % 8;
        assert(byte_index < T::NUM_LETTER_BYTES);
        uint8_t bit = ONE_8;
        bit <<= (7 - bit_left_index);
        data.letters[byte_index] |= bit;
    }

    vec_ind_pair_t get_letter_and_node_indices_for_on_bits() const  {
        assert(!is_terminal());
        vec_ind_pair_t ret;
        ret.reserve(8*T::NUM_LETTER_BYTES - 2);
        u_int32_t node_index = get_index();
        uint8_t masked = data.letters[0] & TOP_LETTER_MASK8;
        fill_letter_and_node_indices(masked, -2, ret, node_index);
        int offset = 6;
        for (unsigned char i = 1; i < T::NUM_LETTER_BYTES; ++i) {
            fill_letter_and_node_indices(data.letters[i], offset, ret, node_index);
            offset += 8;
        }
        return ret;
    }

};

using CTrie128Node = CTrieNode<CTrie128NodeData>;
using CTrie80Node = CTrieNode<CTrie80NodeData>;

} // namespace otc
#endif

        