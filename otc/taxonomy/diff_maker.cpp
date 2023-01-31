#include <boost/algorithm/string/join.hpp>
#include <filesystem>
namespace fs = std::filesystem;
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
                                     OttId kr,
                                     bool read_syn_type_as_src)
    :RichTaxonomy(dir, cf, kr, read_syn_type_as_src) {
    const auto & rich_tax_tree = this->get_tax_tree();
    const auto & rt_data = rich_tax_tree.get_data();
    std::vector<const RTRichTaxNode *> nd_vec;
    nd_vec.reserve(10);
    for (auto& [name, node] : rt_data.name_to_node) {
        nd_vec.clear();
        if (node == nullptr) {
            nd_vec = rt_data.homonym_to_nodes.at(name);
        } else {
            nd_vec.push_back(node);
        }
        for (auto nd : nd_vec ) {
            const auto & nd_data = nd->get_data();
            for (auto syn_ptr : nd_data.junior_synonyms) {
                synonym2node[syn_ptr->name].push_back(nd);
                //std::cerr << "registered " << syn_ptr->name << " syn for "<< nd_data.get_nonuniqname() << '\n';
            }
        }
    }

}


} //namespace otc
