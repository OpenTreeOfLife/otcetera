// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <bitset>
#include <fstream>
#include <regex>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/join.hpp>
namespace fs = boost::filesystem;

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/patching.h"
#include "otc/config_file.h"
#include "otc/util.h"
#include "otc/otc_base_includes.h"
#include "otc/ctrie/context_ctrie_db.h"

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::ofstream;
using std::map;
using std::set;
using std::optional;
using nlohmann::json;
using std::ifstream;
using std::unordered_set;
using std::string_view;
using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;
namespace po = boost::program_options;
using po::variables_map;

namespace otc
{
PatchableTaxonomy load_patchable_taxonomy(const boost::program_options::variables_map& args) {
    string taxonomy_dir = get_taxonomy_dir(args);
    OttId keep_root = -1;
    bitset<32> cleaning_flags = 0;
    return {taxonomy_dir, cleaning_flags, keep_root};
}

PatchableTaxonomy::PatchableTaxonomy(const std::string& dir,
                                     std::bitset<32> cf,
                                     OttId kr)
    :RichTaxonomy(dir, cf, kr) {
}


} //namespace otc
