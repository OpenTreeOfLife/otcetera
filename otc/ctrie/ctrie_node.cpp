#include "otc/ctrie/ctrie_node.h"

using std::optional;

namespace otc {

void CTrieNode::flag_letter(unsigned int letter)
{
    letter_bits |= (1UL<<letter);
}

optional<std::size_t> CTrieNode::child_index_for_letter(stored_index_t letter) const
{
    auto lbits = get_letter_bits();

    // Quit if we don't have a bit for that letter.
    if ((lbits & (1UL<<letter)) == 0) return {};

    // Get bits for letters BEFORE this letter
    uint64_t prev_bits = (lbits<<(63-letter))<<1;

    // The number of letters BEFORE this letter
    int delta = __builtin_popcountl(prev_bits);

    return get_index() + delta;
}

}
