#include "otc/ctrie/ctrie_node.h"

using std::optional;

namespace otc {

void CTrieNode::flag_letter(unsigned int i)
{
    letter_bits |= HIGHEST_BIT >> i;
}

optional<std::size_t> CTrieNode::child_index_for_letter(stored_index_t letter) const
{
    auto lbits = get_letter_bits();

    // We don't have a bit for that letter.
    if ((lbits & (HIGHEST_BIT>>letter)) == 0) return {};

    // remove bits for previous letters
    uint64_t mask = ((-1UL)<<(63-letter))<<1;

    // The number of letters BEFORE this letter
    int delta = __builtin_popcountl(mask & lbits);

    return get_index() + delta;
}

}
