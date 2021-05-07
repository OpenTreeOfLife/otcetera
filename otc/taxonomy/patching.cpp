#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;
#include "otc/taxonomy/patching.h"
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
    const auto & rich_tax_tree = this->get_tax_tree();
    const auto & rt_data = rich_tax_tree.get_data();
    for (auto& [name, node] : rt_data.name_to_node) {
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


std::pair<bool, std::string> PatchableTaxonomy::add_new_taxon(OttId oid,
                                                         OttId parent_id,
                                                         const std::string & name,
                                                         const std::string & rank,
                                                         const std::string & sourceinfo,
                                                         const std::string & uniqname,
                                                         const std::string & flags,
                                                         OttId * homonym_of) {
    using bool_str_t = std::pair<bool, std::string>;
    auto & rich_tax_tree = this->get_tax_tree();
    auto & rt_data = rich_tax_tree.get_data();
    auto nm_nd_it = rt_data.name_to_node.find(name);
    if (nm_nd_it == rt_data.name_to_node.end()) {
        if (homonym_of != nullptr) {
            return bool_str_t{false, "not a homonym"};
        }
    } else if (homonym_of == nullptr) {
        std::string expl = name;
        expl += " is a homonym of " + std::to_string(nm_nd_it->second->get_ott_id());
        return bool_str_t{false, expl};
    }
    vector<string> elements;
    elements.reserve(8);
    elements.push_back(std::to_string(oid));
    elements.push_back(std::to_string(parent_id));
    elements.push_back(name);
    elements.push_back(rank);
    elements.push_back(sourceinfo);
    elements.push_back(uniqname);
    elements.push_back(flags);
    elements.push_back(string());
    string fake_line = boost::algorithm::join(elements, "\t|\t");
    std::cerr << fake_line << std::endl;
    return bool_str_t{false, "not implemented"}; //add_taxon_record(fake_line);
}

void PatchableTaxonomy::write_version_file_contents(std::ostream & out) const {
    out << version << std::endl;
}

void PatchableTaxonomy::write_taxonomy_file_contents(std::ostream & tf) const {
    tf << "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t" << std::endl;
    const auto & tax_tree = get_tax_tree();
    const string sep = "\t|\t";
    for (auto nd : iter_pre_const(tax_tree)) {
        tf << nd->get_ott_id() << sep;
        const auto par = nd->get_parent();
        if (par) {
            tf << par->get_ott_id();
        }
        tf << sep;
        const auto & data = nd->get_data();
        const auto nu_name = data.get_nonuniqname();
        tf << nu_name << sep ;
        tf << data.get_rank() << sep;
        tf << data.source_info << sep;
        const auto & nname = nd->get_name();
        if (nu_name != nname) {
            tf << nname;
        }
        tf << sep;
        tf << flags_to_string(data.flags) << sep;
        tf << '\n';
    }
    // tf << "name2node\n";
    // const auto & rich_tax_tree = this->get_tax_tree();
    // const auto & rt_data = rich_tax_tree.get_data();
    // for (auto & [name, node] : rt_data.name_to_node) {
    //     tf << name << sep << node->get_ott_id() << '\n';
    // }
    // 
    // for (auto& rec: *this) {
    //    tf << rec.id << sep;
    //    if (rec.parent_id > 0) {
    //        tf << rec.parent_id;
    //    }
    //    tf << sep;
    //    tf << rec.name << sep;
    //    tf << rec.rank << sep;
    //    tf << rec.sourceinfo << sep;
    //    if (rec.uniqname != rec.name) {
    //        tf << rec.uniqname;
    //    }
    //    tf << sep;
    //    tf << flags_to_string(rec.flags) << sep;
    //    tf << '\n';
    // }
}

void PatchableTaxonomy::write_synonyms_file_contents(std::ostream & sf) const {
    const string sep = "\t|\t";
    sf << "name\t|\tuid\t|\ttype\t|\tuniqname\t|\tsourceinfo\t|\t" << std::endl;
    for (auto & [syn_name, vec_nd] : synonym2node) {
        for (auto nd_ptr : vec_nd) {
            const auto & jsv = nd_ptr->get_data().junior_synonyms;
            for (auto jsp : jsv) {
                if (jsp->name == syn_name) {
                    sf << syn_name << sep 
                       << nd_ptr->get_ott_id() << sep 
                       << sep // we don't retain the type on parsing
                       << sep // we don't retain the uniqname on parsing
                       << jsp->source_string << sep << '\n';
                    break;
                }
            }
        }
    }
    sf.flush();
}

void PatchableTaxonomy::write_forwards_file_contents(std::ostream & ff) const {
    ff << "id\treplacement\n";
    for(const auto& p: forwards) {
        ff << p.first << '\t' << p.second << '\n';
    }
}


void PatchableTaxonomy::write(const std::string& newdirname) const{
    fs::path old_dir = path;
    fs::path new_dir = newdirname;
    if (! fs::exists(new_dir)) {
        fs::create_directories(new_dir);
    }
    {
        ofstream version_file((new_dir/"version.txt").string());
        write_version_file_contents(version_file);
        version_file.close();
    }
    {
        ofstream tf ((new_dir/"taxonomy.tsv").string());
        write_taxonomy_file_contents(tf);
        tf.close();
    }
    {
        ofstream sf((new_dir/"synonyms.tsv").string());
        write_synonyms_file_contents(sf);
        sf.close();
    }
    {
        ofstream ff((new_dir/"forwards.tsv").string());
        write_forwards_file_contents(ff);
        ff.close();
    }
}

void PatchableTaxonomy::write_to_stream(std::ostream & out) const {
    out << "version.txt\n";
    write_version_file_contents(out);
    out << "taxonomy.tsv\n";
    write_taxonomy_file_contents(out);
    out << "synonyms.tsv\n";
    write_synonyms_file_contents(out);
    out << "forwards.tsv\n";
    write_forwards_file_contents(out);
}

} //namespace otc
