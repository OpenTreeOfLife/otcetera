#ifndef OTC_CTRIE_NODE_H
#define OTC_CTRIE_NODE_H
#include <cstddef>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <bitset>
#include <optional>
#include "otc/error.h"
#include "otc/otc_base_includes.h"

namespace otc {

using stored_index_t = unsigned char;
using ind_pair_t = std::pair<stored_index_t, uint64_t>;

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
constexpr uint64_t ZERO_64 = 0;
constexpr uint64_t ONE_64 = 1;
constexpr uint64_t HIGHEST_BIT = ONE_64 << 63;
constexpr uint64_t SECOND_HIGHEST_BIT = ONE_64 << 62;
constexpr uint64_t INDEX_MASK = ~(HIGHEST_BIT | SECOND_HIGHEST_BIT);
constexpr stored_index_t NO_MATCHING_CHAR_CODE = 255;

// This class is the type of children.end()
// It exists only to be compared against iterators.
class ctrie_child_sentinel
{
};

// This iterator walks forward through the children of a given node.
// The current letter is indicated by the highest non-zero bit in `letter_bits`.
// Bits are numbered in reverse order from normal, with the highest bit (bit 63) indicating letter 0;
class ctrie_child_iterator
{
    uint64_t letter_bits;
    uint64_t index_;

    void mask_cur_letter()
    {
        assert(not done());
        uint64_t curr_bit = 1UL<<letter();
        letter_bits &= (~curr_bit);
    }

    // If there are no 1 bits, then we have visited all the children.
    bool done() const {return not letter_bits;}

public:
    stored_index_t letter() const {assert(not done()); return __builtin_ctzl(letter_bits);}
    uint64_t index() const {assert(not done()); return index_;}

    ctrie_child_iterator operator++() {index_++; mask_cur_letter(); return (*this);}
    ctrie_child_iterator operator++(int) {auto tmp = *this; (*this)++; return tmp;}

    ind_pair_t operator*() const { return ind_pair_t(letter(),index());}

    bool operator==(const ctrie_child_iterator& i2) {return index() == i2.index() and letter() == i2.letter();}
    bool operator!=(const ctrie_child_iterator& i2) {return not (*this == i2);}
    bool operator==(const ctrie_child_sentinel&) {return done();}
    bool operator!=(const ctrie_child_sentinel&) {return not done();}

    ctrie_child_iterator(uint64_t ul, uint64_t ui): letter_bits(ul),index_(ui) {}
};

// This is the range object with begin() and end() methods for use in
// range-for loops: for(auto [letter,index] : nd->children())
struct ctrie_children
{
    ctrie_child_iterator begin_;
    ctrie_child_iterator begin() const {return begin_;}
    ctrie_child_sentinel end() const {return {};}
    ctrie_children(uint64_t ul, uint64_t ui):begin_(ul,ui) {}
};

using vec_ind_pair_t = std::vector<ind_pair_t>;

class CTrieNode {
private:
    uint64_t letter_bits = 0;
    uint64_t index_bits = 0;
    

    void db_write_state(std::ostream &out) const {
       out << "letter_bits=" << std::bitset<64>{letter_bits}
           << " index_bits=" << std::bitset<64>{index_bits} ;
    }

public:

    CTrieNode() { }

    uint64_t get_index() const {
        return INDEX_MASK & index_bits;
    }

    uint64_t get_letter_bits() const
    {
        return letter_bits;
    }

    void log_state() const {
       std::cerr << " CTrieNode( "; db_write_state(std::cerr); std::cerr << ")\n";
    }

    ctrie_children children() const
    {
        return {get_letter_bits(),get_index()};
    }

    void flag_as_key_terminating()
    {
        index_bits |= SECOND_HIGHEST_BIT;
    }

    bool is_key_terminating() const
    {
        return index_bits & SECOND_HIGHEST_BIT;
    }

    void flag_as_terminal()
    {
        index_bits |= HIGHEST_BIT;
    }

    bool is_terminal() const
    {
        return index_bits & HIGHEST_BIT;
    }

    void set_index(std::size_t index)
    {
        index_bits &= ~INDEX_MASK; // sets to 0 any bits for the index
        index_bits |=  (index & INDEX_MASK);
        if (get_index() != index)
            throw OTCError() << "not enough index field to hold pos = " << index;
    }

    void flag_as_suffix(std::size_t pos)
    {
        flag_as_terminal();
        set_index(pos);
    }

    void set_first_child_index(std::size_t index)
    {
        assert(!is_terminal());
        set_index(index);
    }

    std::optional<std::size_t> child_index_for_letter(stored_index_t letter) const;

    void flag_letter(unsigned int i);
};

} // namespace otc
#endif
