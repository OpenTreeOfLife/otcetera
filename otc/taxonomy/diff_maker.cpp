#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;
#include "otc/taxonomy/diff_maker.h"
#include "otc/otc_base_includes.h"

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
// namespace po = boost::program_options;
// using po::variables_map;

namespace otc
{

TaxonomyDiffMaker load_taxonomy_diff_maker(const boost::program_options::variables_map& args) {
    string taxonomy_dir = get_taxonomy_dir(args);
    OttId keep_root = -1;
    bitset<32> cleaning_flags = 0;
    return {taxonomy_dir, cleaning_flags, keep_root};
}

TaxonomyDiffMaker::TaxonomyDiffMaker(const std::string& dir,
                                     std::bitset<32> cf,
                                     OttId kr)
    :RichTaxonomy(dir, cf, kr) {
    const auto & rich_tax_tree = this->get_tax_tree();
    const auto & rt_data = rich_tax_tree.get_data();
    for (auto& [name, node] : rt_data.name_to_node) {
        if (node == nullptr) {
            continue;
        }
        const auto & nd_data = node->get_data();
        for (auto syn_ptr : nd_data.junior_synonyms) {
            synonym2node[syn_ptr->name].push_back(node);
            //std::cerr << "registered " << syn_ptr->name << " syn for "<< nd_data.get_nonuniqname() << '\n';
        }
        //     if (name != node->get_data().get_nonuniqname()) {
        //         assert(!contains(node_to_uniqname, node));
        //         node_to_uniqname[node] = name;
        //     }
    }
    //std::cerr << filtered_records.size() << " filtered_records" << std::endl;
}


} //namespace otc
