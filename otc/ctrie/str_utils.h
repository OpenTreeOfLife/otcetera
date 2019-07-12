#ifndef OTC_CTRIE_STR_UTILS_H
#define OTC_CTRIE_STR_UTILS_H

#include <codecvt>
#include <locale>
#include "otc/error.h"
#include "otc/otc_base_includes.h"

namespace otc {
std::u32string to_u32string(const std::string_view & undecoded);
std::u32string to_u32string(const std::string & undecoded);

#if defined(U32_TRIE_QUERIES)
    using stored_str_t = std::u32string;
    using stored_char_t = char32_t;
    inline stored_str_t to_stored_str_type(const std::string& inp) {
        return to_u32string(inp);
    }
#else
    using stored_str_t = std::string;
    using stored_char_t = char;
    inline stored_str_t to_stored_str_type(const std::string& inp) {
        return inp;
    }
#endif

// you MUST call set_global_conv_facet() before calling any of the functions that use glob_facet
int set_global_conv_facet();

std::string to_char_str(const std::u32string & undecoded);
std::string to_char_str(const char32_t & undecoded);
template<typename T>
void append_char_str(std::string & x, const T & undecoded);
std::string to_char_str(const std::string & undecoded);
std::string to_char_str(const char & undecoded);

// 
extern const std::ctype<char> * glob_facet;


//conversion from https://en.cppreference.com/w/cpp/locale/codecvt
// utility wrapper to adapt locale-bound facets for wstring/wbuffer convert
template <class Facet>
struct deletable_facet : Facet {
    template <class ...Args>
    deletable_facet(Args&& ...args) : Facet(std::forward<Args>(args)...) {}
    ~deletable_facet() {}
};

extern std::wstring_convert<deletable_facet<std::codecvt<char32_t, char, std::mbstate_t> >, char32_t> glob_conv32;
extern std::wstring_convert<std::codecvt_utf8_utf16<char32_t>, char32_t> glob_conv8;

inline std::u32string to_u32string(const std::string_view & undecoded) {
    try {
        return glob_conv32.from_bytes(undecoded.data(), undecoded.data() + undecoded.length());
    } catch (...) {
        throw OTCError() << "Error converting \"" << undecoded << "\" to_u32string";  
    }
}

inline std::u32string to_u32string(const std::string & undecoded) {
    try {
        return glob_conv32.from_bytes(undecoded.data(), undecoded.data() + undecoded.length());
    } catch (...) {
        throw OTCError() << "Error converting \"" << undecoded << "\" to_u32string";  
    }
}


inline std::string to_char_str(const char32_t * undecoded) {
    return glob_conv8.to_bytes(undecoded);
}
 
inline std::string to_char_str(const std::u32string & undecoded) {
    return glob_conv8.to_bytes(undecoded);
}

inline std::string to_char_str(const char32_t & undecoded) {
    return glob_conv8.to_bytes(undecoded);
}


template<>
inline void append_char_str(std::string & x, const char32_t & undecoded) {
    x += to_char_str(undecoded);
}

template<>
inline void append_char_str(std::string & x, const char & undecoded) {
    x.push_back(undecoded);
}

inline std::string to_char_str(const char *  undecoded) {
    return std::string{undecoded};
}

inline std::string to_char_str(const std::string & undecoded) {
    return undecoded;
}

inline std::string to_char_str(const char & undecoded) {
    return std::string{1, undecoded};
}


inline std::u32string to_u32string_ci(const std::string_view & uncap_mod) {
    assert(glob_facet != nullptr);
    std::string undecoded{uncap_mod};
    glob_facet->tolower(&undecoded[0], &undecoded[0] + undecoded.size());
    std::u32string ret;
    try {
        ret = glob_conv32.from_bytes(undecoded.data(), undecoded.data() + undecoded.length());
    } catch (...) {
        throw OTCError() << "Error converting \"" << undecoded << "\" to_u32string";  
    }

    return ret;
}

inline std::string lower_case_version(const std::string & arg) {
    assert(glob_facet != nullptr);
    std::string undecoded{arg};
    glob_facet->tolower(&undecoded[0], &undecoded[0] + undecoded.size());
    return undecoded;
}

inline std::string upper_case_version(const std::string & arg) {
    assert(glob_facet != nullptr);
    std::string undecoded{arg};
    glob_facet->toupper(&undecoded[0], &undecoded[0] + undecoded.size());
    return undecoded;
}


inline bool starts_with(const stored_str_t & full, const stored_str_t & pref) {
    if (full.length() < pref.length()) {
        return false;
    }
    //std::cerr << "starts_with(" << to_char_str(full) << ", " << to_char_str(pref) << ")\n";
    return 0 == full.compare(0, pref.length(), pref);
}

void init_char_maps(); // must be called before normalize_queries are used.
// could use this for \"e  -> e as well
std::string normalize_query(const std::string & raw_query);
std::string normalize_query(const std::string_view & raw_query);



} // namespace otc
#endif
