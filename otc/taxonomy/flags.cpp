#include "otc/taxonomy/flags.h"

// TODO: write out a reduced taxonomy

#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <boost/algorithm/string/join.hpp>
#include <bitset>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/flags.h"

using namespace otc;

using std::string;
using std::vector;
using std::bitset;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

// https://github.com/OpenTreeOfLife/taxomachine/blob/master/src/main/java/org/opentree/taxonomy/OTTFlag.java

// http://www.boost.org/doc/libs/1_50_0/libs/spirit/doc/html/spirit/qi/reference/string/symbols.html

auto get_symbols() {
    symbols<char, int> sym;
    sym.add
        ("not_otu", 0)
        ("environmental", 1)
        ("environmental_inherited", 2)
        ("viral", 3)
        ("hidden", 4)
        ("hidden_inherited", 5)
//        unclassified_direct
        ("was_container", 6)
        ("barren", 8)
        ("extinct", 9)
//        extinct_direct)
        ("extinct_inherited", 11)
        ("major_rank_conflict", 12)
//        major_rank_conflict_direct
        ("major_rank_conflict_inherited", 14)
        ("unclassified", 15)
        ("unclassified_inherited", 16)
        ("edited", 17)
        ("hybrid", 18)
        ("incertae_sedis", 19)
        ("incertae_sedis_inherited", 20)
//     incertae_sedis_direct
        ("infraspecific", 22)
        ("sibling_lower", 23)
        ("sibling_higher", 24)
        ("tattered", 25)
        ("tattered_inherited", 26)
        ("forced_visible", 27)
        ("unplaced", 28)
        ("unplaced_inherited", 29)
        ("inconsistent", 30)
        ("merged", 31)
        ;
    return sym;
}

auto flag_symbols = get_symbols();

namespace otc {
int flag_from_string(const char* start, const char* end) {
    int n = end - start;
    if (n == 0) {
        throw OTCError() << "Flags string with consecutive commas or a starting or trailing commas.";
    }
    assert(n > 0);
    int flag = 0;
    auto cur = start;
    boost::spirit::qi::parse(cur, end, flag_symbols, flag);
    if (cur != end) {
        throw OTCError() << "Flag '" << string(start, end) << "' not recognized.";
    }
    return flag;
}

int flag_from_string(const string& s) {
    const char* start = s.c_str();
    const char* end = start + s.length();
    return flag_from_string(start, end);
}

tax_flags flags_from_string(const char* start, const char* end) {
    assert(start <= end);
    bitset<32> flags;
    while (start < end) {
        assert(start <= end);
        const char* sep = std::strchr(start, ',');
        if (not sep) {
            sep = end;
        }
        int flag = flag_from_string(start, sep);
        flags |= (1 << flag);
        start = sep + 1;
    }
    return flags;
}

tax_flags flags_from_string(const string& s) {
    if (s.empty()) {
        return {};
    }
    const char* start = s.c_str();
    const char* end = start + s.length();
    return flags_from_string(start, end);
}

tax_flags cleaning_flags_from_config_file(const string& filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(filename, pt);
    string cleaning_flags_string = pt.get<std::string>("taxonomy.cleaning_flags");
    return flags_from_string(cleaning_flags_string);
}

tax_flags regrafting_flags_from_config_file(const std::string& filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(filename, pt);
    auto cfs =  pt.get<std::string>("taxonomy.cleaning_flags");
    tax_flags cf;
    if (!cfs.empty()) {
        cf = flags_from_string(cfs);
    }
    tax_flags arf;
    string ars = pt.get<std::string>("taxonomy.additional_regrafting_flags");
    if (!ars.empty()) {
          arf = flags_from_string(ars);
    }
    tax_flags u = cf | arf;
    return u;
}

string string_for_flag(int i) {
    vector<string> matches;
    flag_symbols.for_each([&](const string& s, int j) {
                             if (i == j) {
                                 matches.push_back(s);
                             }
                         });
    return matches[0];
}

std::string flags_to_string(const tax_flags flags) {
    vector<string> f = flags_to_string_vec(flags);
    return boost::algorithm::join(f, ", ");
}

std::vector<std::string> flags_to_string_vec(const std::bitset<32> flags) {
    vector<string> f;
    for(int i=0;i<32;i++) {
        if (flags.test(i)) {
            f.push_back(string_for_flag(i));
        }
    }
    return f;
}

} // namespace otc
