#ifndef OTCETERA_UTIL_H
#define OTCETERA_UTIL_H

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <list>
#include "otc/otc_base_includes.h"
#include "otc/error.h"

namespace otc {

enum QuotingRequirementsEnum {
    NO_QUOTES_NEEDED,
    UNDERSCORE_INSTEAD_OF_QUOTES,
    QUOTES_NEEDED
};

const std::string read_str_content_of_utf8_file(const std::string &filepath);
bool open_utf8_file(const std::string &filepath, std::ifstream & inp);
std::list<std::string> read_lines_of_file(const std::string & filepath);
std::string filepath_to_filename(const std::string &filepath);

bool char_ptr_to_long(const char *c, long *n);
std::size_t find_first_graph_index(const std::string & s);
std::size_t find_last_graph_index(const std::string & s);
std::string strip_leading_whitespace(const std::string & s);
std::string strip_trailing_whitespace(const std::string & s);
std::string strip_surrounding_whitespace(const std::string &n);
// w/o delimiter split_string splits on whitespace and does not return empty 
//  elements. With the delimiter, consecutive delimiters will lead to an empty string
std::list<std::string> split_string(const std::string &s);
std::list<std::string> split_string(const std::string &s, const char delimiter);
std::list<std::set<long> > parse_designators_file(const std::string &fp);
std::set<long> parse_list_of_ott_ids(const std::string &fp);

QuotingRequirementsEnum determine_newick_quoting_requirements(const std::string & s);
std::string add_newick_quotes(const std::string &s);
std::string blanks_to_underscores(const std::string &s);
void write_escaped_for_newick(std::ostream & out, const std::string & n);
void write_ott_id_set(std::ostream & out, const char *indent, const std::set<long> &fir, const char * sep);
void db_write_ott_id_set(const char * label, const std::set<long> &fir);
void write_ott_id_set_diff(std::ostream & out, const char *indent, const std::set<long> &fir, const char *firN, const std::set<long> & sec, const char *secN);
template<typename T>
std::set<T> container_as_set(const std::vector<T> &);
template<typename T>
bool is_subset(const T & small, const T & big);
template<typename T, typename U>
bool have_intersection(const T & first, const U & second);
template<typename T>
std::set<T> set_intersection_as_set(const std::set<T> & small, const std::set<T> & big);
template<typename T>
std::set<T> set_union_as_set(const std::set<T> & small, const std::set<T> & big);
template<typename T>
std::set<T> set_sym_difference_as_set(const std::set<T> & small, const std::set<T> & big);
template<typename T>
std::set<T> set_difference_as_set(const std::set<T> & small, const std::set<T> & big);
template<typename T, typename U>
bool contains(const T & container, const U & key);
template<typename T>
inline bool vcontains(const std::vector<T> & container, const T & key);
template<typename T, typename U>
std::set<T> keys(const std::map<T, U> & container);
std::set<long> parse_delim_separated_ids(const std::string &str, const char delimiter);

void append_include_leaf_set_as_newick(const char *fn, const OttIdSet & inc, const OttIdSet & ls);
template <typename T>
std::ostream& write_separated_collection(std::ostream& o, const std::set<T>& s, const char * sep);
template <typename T>
std::ostream& write_separated_collection(std::ostream& o, const std::list<T>& s, const char * sep);
template <typename T>
std::ostream& write_separated_collection(std::ostream& o, const std::vector<T>& s, const char * sep);


template<typename T>
inline std::string get_contested_preamble_from_name(const T & nd, const std::string & treeName) {
    std::ostringstream ss;
    ss << nd.get_ott_id() << " \"" << nd.get_name() << "\" contested by \"" << treeName << "\"";
    return ss.str();
}

template<typename T, typename U>
inline std::string get_contested_preamble(const T & nd, const U & tree) {
    return get_contested_preamble_from_name(nd, tree.get_name());
}

#if defined WIDE_STR_VERSION
const std::wstring read_wstr_content_of_utf8_file(const std::string &filepath);
bool openUTF8WideFile(const std::string &filepath, std::wifstream & inp);

inline const std::wstring read_wstr_content_of_utf8_file(const std::string &filepath) {
    std::wifstream inp;
    if (!openUTF8WideFile(filepath, inp)) {
        throw OTCError("Could not open \"" + filepath + "\"");
    }
    const std::wstring utf8content((std::istreambuf_iterator<wchar_t>(inp) ), (std::istreambuf_iterator<wchar_t>()));
    return utf8content;
}
#endif //defined WIDE_STR_VERSION

template<typename T, typename U>
inline bool contains(const T & container, const U & key) {
    return container.find(key) != container.end();
}
template<typename T>
inline bool vcontains(const std::vector<T> & container, const T & key) {
    for (const auto & el : container) {
        if (el == key) {
            return true;
        }
    }
    return false;
}
template<typename T, typename U>
inline std::set<T> keys(const std::map<T, U> & container) {
    std::set<T> k;
    for (const auto & x : container) {
        k.insert(x.first);
    }
    return k;
}

inline void write_ott_id_set(std::ostream & out,
                        const char *indent,
                        const std::set<long> &fir,
                        const char * sep) {
    for (auto rIt = fir.begin(); rIt != fir.end(); ++rIt) {
        if (rIt != fir.begin()) {
            out << sep;
        }
        out << indent << "ott" << *rIt;
    }
}
inline void db_write_ott_id_set(const char * label,
                          const std::set<long> &fir) {
    if (!debugging_output_enabled) {
        return;
    }
    std::cerr << label;
    write_ott_id_set(std::cerr, " ", fir, " ");
    std::cerr << std::endl;
}
inline void write_ott_id_set_diff(std::ostream & out,
                            const char *indent,
                            const std::set<long> &fir,
                            const char *firN,
                            const std::set<long> & sec,
                            const char *secN) {
    for (const auto & rIt : fir) {
        if (sec.find(rIt) == sec.end()) {
            out << indent << "ott" << rIt << " is in " << firN << " but not " << secN << "\n";
        }
    }
    for (const auto & rIt : sec) {
        if (fir.find(rIt) == fir.end()) {
            out << indent << "ott" << rIt << " is in " << secN << " but not " << firN << "\n";
        }
    }
}

template<typename T>
inline bool is_proper_subset(const T & small, const T & big) {
    if (big.size() <= small.size()) {
        return false;
    }
    for (const auto & rIt : small) {
        if (big.find(rIt) == big.end()) {
            return false;
        }
    }
    return true;
}

template<typename T>
inline bool is_subset(const T & small, const T & big) {
    if (big.size() < small.size()) {
        return false;
    }
    for (const auto & rIt : small) {
        if (big.find(rIt) == big.end()) {
            return false;
        }
    }
    return true;
}

// http://stackoverflow.com/posts/1964252/revisions
template<class Set1, class Set2>
inline bool are_disjoint(const Set1 & set1, const Set2 & set2) {
    if (set1.empty() || set2.empty()) {
        return true;
    }
    auto it1 = set1.begin();
    const auto it1End = set1.end();
    auto it2 = set2.begin();
    const auto & it2End = set2.end();
    if (*it1 > *set2.rbegin() || *it2 > *set1.rbegin()) {
        return true;
    }
    while (it1 != it1End && it2 != it2End) {
        if (*it1 == *it2) {
            return false;
        }
        if(*it1 < *it2) {
            it1++;
        } else {
            it2++;
        }
    }
    return true;
}

// http://stackoverflow.com/posts/1964252/revisions
template<typename T, typename U>
inline bool ds_are_disjoint(const std::map<T, U> & first,
                          const std::set<T> & second) {
    if (first.empty() || second.empty()) {
        return true;
    }
    auto it1 = first.begin();
    const auto it1End = first.end();
    auto it2 = second.begin();
    const auto & it2End = second.end();
    if (it1->first > *second.rbegin() || *it2 > (first.rbegin()->first)) {
        return true;
    }
    while (it1 != it1End && it2 != it2End) {
        if (it1->first == *it2) {
            return false;
        }
        if(it1->first < *it2) {
            it1++;
        } else {
            it2++;
        }
    }
    return true;
}


// called when we've determined that set1 is smaller than set2, and they have an
// intersection with the first el of set1, so set 1 could be a subset of set2
// returns true if set1 is a subset of set2 (where the args are the iterators and
// end iterators for each set)
template<typename T>
inline bool finish_subset_compat(T & it1, const T & it1End, T & it2, const T &it2End) {
    while (it1 != it1End && it2 != it2End) {
        if (*it1 == *it2) {
            ++it1;
            ++it2;
        } else if (*it1 < *it2) {
            return false;
        } else {
            ++it2;
        }
    }
    return it1 == it1End;
}

// adapted http://stackoverflow.com/posts/1964252/revisions
template<typename T>
inline bool are_compatible_des_id_sets(const T & set1, const T & set2) {
    if (set1.size() < 2 || set2.size() < 2) {
        return true;
    }
    auto it1 = set1.begin();
    const auto it1End = set1.end();
    auto it2 = set2.begin();
    const auto it2End = set2.end();
    if (*it1 > *set2.rbegin() || *it2 > *set1.rbegin()) {
        return true; // disjoint
    }
    while (it1 != it1End && it2 != it2End) {
        if (*it1 == *it2) {
            // there is an intersection. they must be equal, or one a subset of the other...
            if (it1 == set1.begin()) {
                if (it2 == set2.begin()) {
                    //if they are the same size, they have to be equal
                    if (set1.size() == set2.size()) {
                        return set1 == set2;
                    }
                    ++it1;
                    ++it2;
                    // otherwise, they are only compatible if the
                    //  smaller set is a subset of the other.
                    if (set1.size() < set2.size()) {
                        return finish_subset_compat(it1, it1End, it2, it2End);
                    } else {
                        return finish_subset_compat(it2, it2End, it1, it1End);
                    }
                } else {
                    // did not match first el of set2, so only compat if set1 is in set2
                    if (set2.size() < set1.size()) {
                        return false;
                    }
                    return finish_subset_compat(it1, it1End, it2, it2End);
                }
            } else if (it2 == set2.begin()) {
                // did not match first el of set1, so only compat if set2 is in set1
                if (set1.size() < set2.size()) {
                    return false;
                }
                return finish_subset_compat(it2, it2End, it1, it1End);
            } else {
                // first intersection is not the first element of either,
                //   so one can't be a subset of the other...
                return false;
            }
        }
        if(*it1 < *it2) {
            it1++;
        } else {
            it2++;
        }
    }
    return true; // disjoint, so compatible
}


template<typename T>
inline bool have_intersection(const T & first, const T & second) {
    return !are_disjoint<T>(first, second);
}

template<typename T, typename U>
inline bool ds_have_intersection(const std::map<T, U> & first,
                               const std::set<T> & second) {
    return !ds_are_disjoint<T>(first, second);
}

template<typename T>
inline std::set<T> container_as_set(const std::vector<T> &v) {
    std::set<T> d;
    d.insert(v.begin(), v.end());
    return d;
}
template<typename T>
inline std::set<T> set_intersection_as_set(const std::set<T> & fir, const std::set<T> & sec) {
    std::set<T> d;
    set_intersection(begin(fir), end(fir), begin(sec), end(sec), std::inserter(d, d.end()));
    return d;
}
template<typename T>
inline std::set<T> set_union_as_set(const std::set<T> & fir, const std::set<T> & sec) {
    std::set<T> d;
    set_union(begin(fir), end(fir), begin(sec), end(sec), std::inserter(d, d.end()));
    return d;
}
template<typename T>
inline std::set<T> set_sym_difference_as_set(const std::set<T> & fir, const std::set<T> & sec) {
    std::set<T> d;
    set_symmetric_difference(begin(fir), end(fir), begin(sec), end(sec), std::inserter(d, d.end()));
    return d;
}
template<typename T>
inline std::set<T> set_difference_as_set(const std::set<T> & fir, const std::set<T> & sec) {
    std::set<T> d;
    set_difference(begin(fir), end(fir), begin(sec), end(sec), std::inserter(d, d.end()));
    return d;
}

inline std::size_t find_first_graph_index(const std::string & s) {
    std::size_t pos = 0U;
    for (const auto & c : s) {
        if (isgraph(c)) {
            return pos;
        }
        ++pos;
    }
    return std::string::npos;
}

inline std::size_t find_last_graph_index(const std::string & s) {
    auto pos = s.length();
    while (pos > 0) {
        --pos;
        if (isgraph(s[pos])) {
            return pos;
        }
    }
    return std::string::npos;
}

inline std::string strip_leading_whitespace(const std::string & n) {
    const auto x = find_first_graph_index(n);
    if (x == std::string::npos) {
        return std::string();
    }
    return n.substr(x);
}

inline std::string strip_trailing_whitespace(const std::string & n) {
    const auto x = find_last_graph_index(n);
    if (x == std::string::npos) {
        return std::string();
    }
    return n.substr(0, 1 + x);
}

inline std::string strip_surrounding_whitespace(const std::string &n) {
    const auto s = find_first_graph_index(n);
    if (s == std::string::npos) {
        return std::string();
    }
    const auto e = find_last_graph_index(n);
    assert(e != std::string::npos);
    return n.substr(s, 1 + e - s);
}


inline QuotingRequirementsEnum determine_newick_quoting_requirements(const std::string & s) {
    QuotingRequirementsEnum nrq = NO_QUOTES_NEEDED;
    for (const auto & c : s) {
        if (!isgraph(c)) {
            if (c != ' ') {
                return QUOTES_NEEDED;
            }
            nrq  = UNDERSCORE_INSTEAD_OF_QUOTES;
        } else if (strchr("(){}\"-]/\\,;:=*`+<>", c) != nullptr) {
            return (s.length() > 1 ? QUOTES_NEEDED : NO_QUOTES_NEEDED);
        } else if (strchr("\'[_", c) != nullptr) {
            return QUOTES_NEEDED;
        }
    }
    return nrq;
}

inline std::string add_newick_quotes(const std::string &s) {
    std::string withQuotes;
    unsigned len = static_cast<unsigned>(s.length());
    withQuotes.reserve(len + 4);
    withQuotes.append(1,'\'');
    for (const auto & c : s) {
        withQuotes.append(1, c);
        if (c == '\'') {
            withQuotes.append(1,'\'');
        }
    }
    withQuotes.append(1,'\'');
    return withQuotes;
}

inline std::string blanks_to_underscores(const std::string &s) {
    std::string r{s};
    std::replace(begin(r), end(r), ' ', '_');
    return r;
}

inline void write_escaped_for_newick(std::ostream & out, const std::string & n) {
    const QuotingRequirementsEnum r = determine_newick_quoting_requirements(n);
    if (r == NO_QUOTES_NEEDED) {
        out << n;
    } else if (r == UNDERSCORE_INSTEAD_OF_QUOTES) {
        out << blanks_to_underscores(n);
    } else {
        out << add_newick_quotes(n);
    }
}


template<typename T>
inline T intersection_of_sets(const T & first, const T &sec) {
    T intersection;
    std::set_intersection(first.begin(), first.end(),
                          sec.begin(), sec.end(),
                          std::inserter(intersection, intersection.begin()));
    return intersection;
}
template<typename T>
inline std::size_t size_of_symmetric_difference(const T & first, const T &sec) {
    T diff;
    std::set_symmetric_difference(first.begin(), first.end(),
                                  sec.begin(), sec.end(),
                                  std::inserter(diff, diff.begin()));
    return diff.size();
}

template <typename T>
inline std::ostream& write_separated_collection(std::ostream& o, const std::set<T>& s, const char * sep) {
    auto it = s.begin();
    o << *it++;
    for(; it != s.end(); it++) {
        o << sep << *it;
    }
    return o;
}

template <typename T>
inline std::ostream& write_separated_collection(std::ostream& o, const std::list<T>& s, const char * sep) {
    auto it = s.begin();
    o << *it++;
    for(; it != s.end(); it++) {
        o << sep << *it;
    }
    return o;
}

template <typename T>
inline std::ostream& write_separated_collection(std::ostream& o, const std::vector<T>& s, const char * sep) {
    auto it = s.begin();
    o << *it++;
    for(; it != s.end(); it++) {
        o << sep << *it;
    }
    return o;
}


} //namespace otc
#endif
