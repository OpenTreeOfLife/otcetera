#include "otc/ctrie/ctrie_node.h"

namespace otc {

//@TODO use lookup table
void fill_letter_and_node_indices(unsigned char curr_byte,
                                  int offset,
                                  vec_ind_pair_t & ret,
                                  uint64_t & node_index)
{
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

void fill_letter_and_node_indices_64(uint64_t masked_workspace,
                                     int offset,
                                     vec_ind_pair_t & ret,
                                     uint64_t & node_index)
{
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

template <>
void CTrieNode<CTrie3NodeData>::flag_letter(unsigned int i) {
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
vec_ind_pair_t CTrieNode<CTrie3NodeData>::get_letter_and_node_indices_for_on_bits() const {
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
void CTrieNode<CTrie2NodeData>::flag_letter(unsigned int i) {
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
vec_ind_pair_t CTrieNode<CTrie2NodeData>::get_letter_and_node_indices_for_on_bits() const {
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
void CTrieNode<CTrie1NodeData>::flag_letter(unsigned int i) {
    uint64_t bit = ONE_64;
    assert(i < DATA_TYPE::END_LETTER_INDEX);
    bit <<= (LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD - 1 - i);
    data.top |= bit;
} 

template <>
vec_ind_pair_t CTrieNode<CTrie1NodeData>::get_letter_and_node_indices_for_on_bits() const {
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

}
