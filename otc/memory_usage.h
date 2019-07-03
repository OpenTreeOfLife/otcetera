#ifndef OTCETERA_MEMORY_USAGE_H
#define OTCETERA_MEMORY_USAGE_H
// utility functions for approximating the memory used by data structures.

#include "otc/otc_base_includes.h"
#include "otc/tree_iter.h"
#include "otc/error.h"
#include "otc/util.h"
#include "otc/debug.h"
#include <vector>
#include <unordered_map>
#include <set>
#include <string>
namespace otc {

using MemoryBookkeeper = std::unordered_map<std::string, std::size_t>;
template<typename T>
std::size_t calc_memory_used(const T &, MemoryBookkeeper &);

template<typename T>
inline std::size_t calc_memory_used(const T &, MemoryBookkeeper &) {
    return sizeof(T);
}


template<>
inline std::size_t calc_memory_used(const OttIdSet &d,
                                    MemoryBookkeeper &) {
    return sizeof(char *) * 2 + sizeof(OttId) * d.size();
}

template<typename T>
inline void write_memory_bookkeeping(T & outstream,
                                     MemoryBookkeeper & mb,
                                     const std::string & tag,
                                     std::size_t total) {
    std::string ts = std::to_string(total);
    std::size_t num_digits = std::max<std::size_t>(10, ts.length());
    outstream << "Memory report for " << tag << ":\n";
    outstream << std::setw (num_digits) << total << " " << "Total\n";
    std::set<std::string> ks;
    for (auto el : mb) {
        ks.insert(el.first);
    }
    for (auto k : ks) {
        outstream << std::setw (num_digits) <<  mb.at(k) << " " << k << "\n";
    }
}

template<>
inline std::size_t calc_memory_used(const std::string &s, MemoryBookkeeper &) {
    //ignores reference counting, but adds an int32_t for a ref count
    std::size_t total = s.capacity() * sizeof(char);
    total += sizeof(char *); // start pointer
    total += sizeof(std::size_t); // capacity
    total += sizeof(int32_t); // approximate - may be STL impl dependent
    return total;
}

template<typename T, typename U>
inline std::size_t calc_memory_used_owned_pair(const std::pair<T*, U*> &p, MemoryBookkeeper &mb) {
    std::size_t total = 2 * sizeof(T *);
    if (p.first) {
        total += calc_memory_used(*p.first, mb);
    }
    if (p.second) {
        total += calc_memory_used(*p.second, mb);
    }
    return total;
}

template<typename E>
inline std::size_t calc_memory_used_by_vector_eqsize(const std::vector<E> & v,
                                                     std::size_t el_size,
                                                     MemoryBookkeeper &) {
    std::size_t total = 0;
    total += sizeof(E *); // start pointer
    total += sizeof(std::size_t); // capacity
    total += v.capacity() * el_size;
    return total;
}

template<typename K, typename V>
inline std::size_t calc_memory_used_by_map_eqsize(const std::unordered_map<K, V> & v, std::size_t el_size, MemoryBookkeeper &) {
    std::size_t total = 0;
    const int num_unused_buckets = v.bucket_count() - v.size();
    if (num_unused_buckets > 0) {
        total += num_unused_buckets * (sizeof(K) + sizeof(V));
    }
    total += v.size() * el_size;
    return total;
}

template<typename K, typename V>
inline std::size_t calc_memory_used_by_map_simple(const std::unordered_map<K, V> & v, MemoryBookkeeper &mb) {
    std::size_t total = 0;
    for (auto el : v) {
        total += calc_memory_used(el.first, mb);
        total += calc_memory_used(el.second, mb);
    }
    return total;
}
template<typename K, typename V>
inline std::size_t calc_memory_used_by_map_simple(const std::map<K, V> & v, MemoryBookkeeper &mb) {
    std::size_t total = 0;
    for (auto el : v) {
        total += calc_memory_used(el.first, mb);
        total += calc_memory_used(el.second, mb);
    }
    return total;
}


template<typename N>
inline std::size_t calc_memory_used_by_node(const N & node, MemoryBookkeeper &mb) {
    std::size_t total = 0;
    total += 5 * sizeof(N *); // 2 child, 2 sib and par pointers
    mb["node navigation"] += total;
    std::size_t nn = calc_memory_used(node.get_name(), mb);
    mb["node name"] += nn;
    std::size_t nd = calc_memory_used(node.get_data(), mb);
    mb["node data"] += nd;
    total += nd + nn;
    total += sizeof(OttId);
    mb["node ids"] += sizeof(OttId);
    return total;
}

template<typename T>
inline std::size_t calc_memory_used_by_tree(const T & tree, MemoryBookkeeper &mb) {
    std::size_t td = calc_memory_used(tree.get_data(), mb);
    std::size_t tn = calc_memory_used(tree.get_name(), mb);
    auto dns = tree.get_detached();
    // approximate
    std::size_t dd = dns.size()*sizeof(int *); // count root pointer
    for (auto nd : dns) {
        assert(nd);
        dd += calc_memory_used_by_node(*nd, mb);
    }
    mb["tree data"] += td;
    mb["tree name"] += tn;
    mb["tree detached nodes"] += dd;
    std::size_t total = dd + td + tn;
    total += sizeof(int *); // count root pointer
    for (auto nd : iter_node_const(tree)) {
        assert(nd);
        total += calc_memory_used_by_node(*nd, mb);
    }
    return total;
}

}// namespace otc
#endif
