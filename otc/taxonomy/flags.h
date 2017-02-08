#ifndef OTC_TAXONOMY_FLAGS_H
#define OTC_TAXONOMY_FLAGS_H

#include <bitset>
#include <string>

namespace otc
{
    typedef std::bitset<32> tax_flags;
    int flag_from_string(const char* start, const char* end);
    int flag_from_string(const std::string&);
    tax_flags flags_from_string(const char* start, const char* end);
    tax_flags flags_from_string(const std::string&);
    std::string flag_to_string(int);
    std::string flags_to_string(const tax_flags);
    tax_flags cleaning_flags_from_config_file(const std::string& filename);
}

#endif
