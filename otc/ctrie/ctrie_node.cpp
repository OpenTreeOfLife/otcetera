#include "otc/ctrie/ctrie_node.h"

namespace otc {

void fill_letter_and_node_indices_64(uint64_t letter_bits,
                                     int offset,
                                     vec_ind_pair_t & ret,
                                     uint64_t & node_index)
{
    constexpr uint64_t left_bit = (ONE_64<<63);
    while (letter_bits)
    {
        int i = __builtin_clzl(letter_bits);

        uint64_t curr_bit = left_bit >> i;
        if (letter_bits & curr_bit)
            ret.push_back(ind_pair_t{offset+i, node_index++});

        letter_bits &= (~curr_bit);
    }
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

    vec_ind_pair_t ret2;


    return ret;
}

}
