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

constexpr std::size_t num_index_bits = 50;

// top:
//      highest bit: is terminal node
//      If that is 0:
//          second highest bit: has key that terminates with this node
//          bits 0 - 61 letter codes

constexpr uint64_t ZERO_64 = 0;
constexpr uint64_t ONE_64 = 1;
constexpr uint64_t HIGHEST_BIT = ONE_64 << 63;
constexpr uint64_t SECOND_HIGHEST_BIT = ONE_64 << 62;
constexpr uint64_t INDEX_MASK = (ONE_64 << num_index_bits) - 1;
constexpr stored_index_t NO_MATCHING_CHAR_CODE = 255;

//std::size_t max_node_index = 0;

template<typename T> void ctrien_set_first_child_index(T& node, std::size_t index);
template<typename T> void ctrien_flag_as_key_terminating(T & node);
template<typename T> void ctrien_flag_as_terminal(T & node);
template<typename T> void ctrien_flag_as_suffix(T & node, std::size_t pos);
template<typename T> void ctrien_ctrien_set_index(T& node, std::size_t index);
template<typename T> bool ctrien_is_key_terminating(const T & node);
template<typename T> bool ctrien_is_terminal(const T & node);
template<typename T> uint64_t ctrien_get_index(const T & node);

template<typename T>
inline void ctrien_flag_as_key_terminating(T & node) {
    node.top |= SECOND_HIGHEST_BIT;
}

template<typename T>
inline bool ctrien_is_key_terminating(const T & node) {
    return node.top & SECOND_HIGHEST_BIT;
}

template<typename T>
inline void ctrien_flag_as_terminal(T & node) {
    node.top |= HIGHEST_BIT;
}

template<typename T>
inline bool ctrien_is_terminal(const T & node) {
    return node.top & HIGHEST_BIT;
}


template<typename T>
uint64_t ctrien_get_index(const T & node) {
    uint64_t masked = node.bot & INDEX_MASK;
    return masked;
}


template<typename T>
inline void ctrien_set_index(T& node, std::size_t index) {
    uint64_t ind = index;
    if ((ind & INDEX_MASK) != ind) {
        throw OTCError() << "not enough index field to hold pos = " << index;
    }
    // if (index > max_node_index) {
    //     max_node_index = index;
    // }
    node.bot = ind;
}

template<typename T>
inline void ctrien_flag_as_suffix(T & node, std::size_t pos) {
    ctrien_flag_as_terminal(node);
    ctrien_set_index(node, pos);
}

template<typename T>
inline void ctrien_set_first_child_index(T& node, std::size_t index) {
    assert((node.top & HIGHEST_BIT) == 0);
    ctrien_set_index(node, index);
}

using ind_pair_t = std::pair<stored_index_t, uint64_t>;
using vec_ind_pair_t = std::vector<ind_pair_t>;

inline void fill_letter_and_node_indices(unsigned char curr_byte,
                                         int offset,
                                         vec_ind_pair_t & ret,
                                         uint64_t & node_index) {
    if (curr_byte == 0) {
        return;
    }                             
    unsigned char curr_bit = 1 << 7;
    for (unsigned char i = 0; i < 8; ++i) {
        //std::cerr << "fill_letter_and_node_indices byte=" << std::hex << (unsigned int)curr_byte << " bit=" << std::hex << (unsigned int)curr_bit << '\n' << std::dec;
        if (curr_byte & curr_bit) {
            ret.push_back(ind_pair_t{i + offset, node_index++});
        }
        curr_bit >>= 1;
    }
}


constexpr unsigned char top_first_byte = 0x03F; // 64 and 128 bit used for top 2 flags
constexpr unsigned char full_byte = 0x0FF;
constexpr unsigned char bot_last_byte = 0x0FC; // bits for 1 and 2 used for upper bits of index;

class CTrie3Node {
    public:
    uint64_t top, mid, bot;
    CTrie3Node() :top{ZERO_64}, mid{ZERO_64}, bot{ZERO_64} {
        //log_state();
    }
    void log_state() const {
       std::cerr << " CTrie2Node( ";
       db_write_state(std::cerr);
       std::cerr << ")\n";
    }
    void db_write_state(std::ostream &out) const {
       out << "top=" << std::bitset<64>{top}
           << " mid=" << std::bitset<64>{mid} 
           << " bot=" << std::bitset<64>{bot} ;
    }

    void flag_letter(unsigned int i) {
        uint64_t bit = 1;
        //log_state();
        if (i < 62) {
            const uint64_t shifted = (bit << (61 - i));
            top |= shifted;
        } else if (i < 126) {
            bit <<= (125 - i);
            mid |= bit;
        } else {
            assert(i < 140);
            bit <<= (139 - i);
            bot |= bit;
        }
    }

    vec_ind_pair_t get_letter_and_node_indices_for_on_bits() const {
        assert(!ctrien_is_terminal(*this));
        vec_ind_pair_t ret;
        ret.reserve(62 + 64 + 14);
        uint64_t node_index = ctrien_get_index(*this);
        int offset = -2;
        int bitshift = 56;
        unsigned char curr_byte = (top >> bitshift) & top_first_byte;
        fill_letter_and_node_indices(curr_byte, offset, ret, node_index);
        for (int i = 0; i < 7; ++i) {
            bitshift -= 8;
            offset += 8;
            curr_byte = (top >> bitshift) & full_byte;
            fill_letter_and_node_indices(curr_byte, offset, ret, node_index);
        }
        // 1 and part of a byte in "bot"
        bitshift = 56;
        for (int i = 0; i < 8; ++i) {
            curr_byte = (mid >> bitshift) & full_byte;
            bitshift -= 8;
            offset += 8;
            fill_letter_and_node_indices(curr_byte, offset, ret, node_index);
        }
        // 1 and part of a byte in "bot"
        bitshift = 56;
        curr_byte = (bot >> bitshift) & full_byte;
        offset += 8;
        fill_letter_and_node_indices(curr_byte, offset, ret, node_index);
        bitshift -= 8;
        curr_byte = (bot >> bitshift) & bot_last_byte;
        offset += 8;
        fill_letter_and_node_indices(curr_byte, offset, ret, node_index);
        return ret;
    }

    const static std::size_t max_num_letters = 3*64 - num_index_bits;
};


class CTrie2Node {
    public:
    uint64_t top, bot;
    CTrie2Node() :top(0),  bot(0) {
    }
    void log_state() const {
       std::cerr << " CTrie2Node( ";
       db_write_state(std::cerr);
       std::cerr << ")\n";
    }
    void db_write_state(std::ostream &out) const {
        out << "top=" << std::bitset<64>{top} 
            << " bot=" << std::bitset<64>{bot};
    }
    
    void flag_letter(unsigned int i) {
        uint64_t bit = 1;
        if (i < 62) {
            bit <<= (61 - i);
            top |= bit;
        }  else {
            assert(i < 76);
            bit <<= (75 - i);
            bot |= bit;
        }
    } 

    vec_ind_pair_t get_letter_and_node_indices_for_on_bits() const {
        //std::cerr << "get_letter_and_node_indices_for_on_bits top="
        //          << std::hex << top << " bot=" << std::hex << bot << std::dec << '\n';
        assert(!ctrien_is_terminal(*this));
        vec_ind_pair_t ret;
        ret.reserve(62 + 64 + 14);
        u_int64_t node_index = ctrien_get_index(*this);
        unsigned char curr_byte = (top >> 56) & top_first_byte;
        int offset = -2;
        fill_letter_and_node_indices(curr_byte, offset, ret, node_index);
        int bitshift = 48;
        for (int i = 0; i < 7; ++i) {
            curr_byte = (top >> bitshift) & full_byte;
            offset += 8;
            fill_letter_and_node_indices(curr_byte, offset, ret, node_index);
            bitshift -= 8;
        }
        // 1 and part of a byte in "bot"
        bitshift = 56;
        curr_byte = (bot >> bitshift) & full_byte;
        offset += 8;
        fill_letter_and_node_indices(curr_byte, offset, ret, node_index);
        bitshift -= 8;
        curr_byte = (bot >> bitshift) & bot_last_byte;
        offset += 8;
        fill_letter_and_node_indices(curr_byte, offset, ret, node_index);
        return ret;
    }

    const static std::size_t max_num_letters = 2*64 - num_index_bits; 
};

} // namespace otc
#endif
