#include "otc/ctrie/ctrie_node.h"

namespace otc {

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

}
