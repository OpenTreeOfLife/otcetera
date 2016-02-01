#ifndef OTC_TAXONOMY_FLAGS_H
#define OTC_TAXONOMY_FLAGS_H

#include <bitset>
#include <string>

namespace otc
{
    std::bitset<32> flag_from_string(const char* start, const char* end);
    std::bitset<32> flag_from_string(const std::string&);
    std::bitset<32> flags_from_string(const char* start, const char* end);
    std::bitset<32> flags_from_string(const std::string&);
    std::bitset<32> cleaning_flags_from_config_file(const std::string& filename);
}

#endif
