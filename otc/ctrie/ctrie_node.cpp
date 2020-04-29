#include "otc/ctrie/ctrie_node.h"

namespace otc {

void CTrieNode::flag_letter(unsigned int i)
{
    uint64_t bit = ONE_64;
    if (i < LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD) {
        bit <<= (LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD - 1 - i);
        top |= bit;
    }  else {
        assert(i < END_LETTER_INDEX);
        bit <<= (LETTER_INDEX_OF_FIRST_BIT_IN_SECOND_WORD - 1 - i);
        bot |= bit;
    }
}

}
