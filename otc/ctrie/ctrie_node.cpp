#include "otc/ctrie/ctrie_node.h"

using std::optional;

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

optional<std::size_t> CTrieNode::child_index_for_letter(stored_index_t letter) const
{
    auto lbits = get_letter_bits();

    // We don't have a bit for that letter.
    if ((lbits & (ONE_64<<(63-letter))) == 0) return {};

    // remove bits for previous letters
    uint64_t mask = ((uint64_t)(-1))<<(64-letter);

    // The number of letters BEFORE this letter
    int delta = __builtin_popcountl(mask & lbits);

    return get_index() + delta;
}

}
