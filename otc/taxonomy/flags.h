#ifndef OTC_TAXONOMY_FLAGS_H
#define OTC_TAXONOMY_FLAGS_H

#include <bitset>
#include <string>
#include <vector>
#include <optional>

namespace otc {

typedef std::bitset<32> tax_flags;
std::optional<std::string> string_for_flag(int i);
int flag_from_string(const char* start, const char* end);
int flag_from_string(const std::string&);
tax_flags flags_from_string(const char* start, const char* end);
tax_flags flags_from_string(const std::string&);
std::string flag_to_string(int);
std::string flags_to_string(const tax_flags);
std::vector<std::string> flags_to_string_vec(const tax_flags flags);
std::bitset<32> cleaning_flags_from_config_file(const std::string& filename);
std::bitset<32> regrafting_flags_from_config_file(const std::string& filename);

bool is_extinct(tax_flags);

} //namespace otc

#endif
