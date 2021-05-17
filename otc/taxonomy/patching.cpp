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

using bool_str_t = std::pair<bool, std::string>;



bool_str_t PatchableTaxonomy::add_synonym(const std::string & name,
                                          OttId ott_id,
                                          const std::string & sourceinfo) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    auto target_nd = included_taxon_from_id(ott_id); 
    if (target_nd == nullptr) {
        auto itrit = rt_data.id_to_record.find(ott_id);
        if (itrit == rt_data.id_to_record.end()) {
            std::string expl = "OTT ID " + std::to_string(ott_id) + " unrecognized.";
            return bool_str_t{false, expl};
        }
        LightSynonym ls{name, sourceinfo};
        auto tr = itrit->second;
        rec_to_new_syn[tr].push_back(ls);
        return bool_str_t{false, ""};
    }
    synonyms.emplace_back(name, target_nd, sourceinfo);
    const TaxonomicJuniorSynonym * syn_ptr = &(*synonyms.rbegin());
    RTRichTaxNodeData & nd_data = const_cast<RTRichTaxNodeData &>(target_nd->get_data());
    nd_data.junior_synonyms.push_back(syn_ptr);
    add_name_to_node_maps(name, target_nd);
    synonym2node[name].push_back(target_nd);
    return bool_str_t{true, ""};
}

template<typename T>
std::vector<T> copy_except(const std::vector<T> & src, const T & taboo) {
    std::vector<T> ret;
    ret.reserve(src.size());
    for (auto i : src) {
        if (i == taboo) {
            continue;
        }
        ret.push_back(i);
    }
    return ret;
}

bool_str_t PatchableTaxonomy::delete_synonym(const std::string & name, OttId ott_id) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    auto target_nd = included_taxon_from_id(ott_id); 
    if (target_nd == nullptr) {
        std::string expl = "OTT ID " + std::to_string(ott_id) + " unrecognized.";
        return bool_str_t{false, expl};
    }
    RTRichTaxNodeData & nd_data = const_cast<RTRichTaxNodeData &>(target_nd->get_data());
    unsigned njsm = 0;
    std::vector<const TaxonomicJuniorSynonym *> njsv;
    for (auto tjs : nd_data.junior_synonyms) {
        if (tjs->name == name) {
            ++njsm;
        } else {
            njsv.push_back(tjs);
        }
    }
    if (njsm == 0) {
        std::string expl = "No synonym of " + std::to_string(ott_id) + " had the name \"" + name + "\".";
        return bool_str_t{false, expl};
    }
    nd_data.junior_synonyms == njsv;
    remove_name_to_node_from_maps(name, target_nd);
    auto & csv = synonym2node[name];
    auto nv = copy_except(csv, target_nd);
    if (nv.size() != csv.size()) {
        synonym2node[name] = nv;
    }
    return bool_str_t{true, ""};

}

void PatchableTaxonomy::add_name_to_node_maps(const std::string & name,
                                                      const RTRichTaxNode * target_nd) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    auto hit = rt_data.homonym_to_nodes.find(name);
    auto sit = rt_data.name_to_node.find(name);
    if (hit != rt_data.homonym_to_nodes.end()) {
        assert(sit == rt_data.name_to_node.end());
        auto nv = copy_except(hit->second, target_nd);
        auto nvs = nv.size();
        if (nvs == hit->second.size()) {
            // if the same size, then target_nd was not in the vector, so add it
            hit->second.push_back(target_nd);
        }
        return;
    }
    if (sit == rt_data.name_to_node.end()) {
        rt_data.name_to_node[name] = target_nd;
        return;
    }
    if (sit->second != target_nd) {
        std::vector<const RTRichTaxNode *> v{sit->second, target_nd};
        rt_data.homonym_to_nodes.emplace(name, v);
        rt_data.name_to_node.erase(sit);
        return;
    }
}


void PatchableTaxonomy::remove_name_to_node_from_maps(const std::string & name,
                                                      const RTRichTaxNode * target_nd) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    auto hit = rt_data.homonym_to_nodes.find(name);
    auto sit = rt_data.name_to_node.find(name);
    if (hit != rt_data.homonym_to_nodes.end()) {
        assert(sit == rt_data.name_to_node.end());
        auto nv = copy_except(hit->second, target_nd);
        auto nvs = nv.size();
        if (nvs == hit->second.size() || nvs == 0) {
            return;
        } else if (nvs == 1) {
            auto nd_p = nv[0];
            rt_data.name_to_node[name] = nd_p;
            return;
        } else {
            rt_data.homonym_to_nodes[name] = nv;
        }
        return;
    }
    if (sit == rt_data.name_to_node.end()) {
        // name must map to taxonrecord instead of a node
        return ;
    }
    if (sit->second != target_nd) {
        return ; // odd
    }
    rt_data.name_to_node.erase(sit);
}

bool_str_t PatchableTaxonomy::add_forward(OttId former_id, OttId redirect_to_id) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    auto fnd = included_taxon_from_id(former_id);
    auto itrit = rt_data.id_to_record.find(former_id);
    if (fnd != nullptr || itrit != rt_data.id_to_record.end()) {
        std::string expl = "OTT ID " + std::to_string(former_id);
        expl += " is in use. It must be deleted before a forward from it is added.";
        return bool_str_t{false, expl};
    }
    fnd = included_taxon_from_id(redirect_to_id);
    itrit = rt_data.id_to_record.find(redirect_to_id);
    if (fnd == nullptr && itrit == rt_data.id_to_record.end()) {
        std::string expl = "OTT ID " + std::to_string(redirect_to_id);
        expl += " not recognized. It cannot be a redirection target.";
        return bool_str_t{false, expl};
    }
    forwards[former_id] = redirect_to_id;
    return bool_str_t{true, ""};
}

bool_str_t PatchableTaxonomy::delete_forward(OttId former_id, OttId redirect_to_id) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    throw OTCError() << "delete_forward not implemented";
}

bool_str_t PatchableTaxonomy::delete_taxon(OttId ott_id) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    throw OTCError() << "delete_taxon not implemented";
}

bool_str_t PatchableTaxonomy::edit_taxon(OttId oid,
                                         OttId parent_id,
                                         const std::string & name,
                                         const std::string & rank,
                                         const std::string & sourceinfo,
                                         const std::string & uniqname,
                                         const std::string & flags,
                                         OttId * homonym_of) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    auto nd_ptr = const_cast<RTRichTaxNode * >(included_taxon_from_id(oid));
    if (nd_ptr == nullptr) {
        std::string expl = "OTT ID " + std::to_string(oid) + " is not an included taxon.";
        return bool_str_t{false, expl};
    }
    RTRichTaxNodeData & nd_data = nd_ptr->get_data();
    if (parent_id != UINT_MAX && parent_id != 0) {
        auto old_par = const_cast<RTRichTaxNode * >(included_taxon_from_id(parent_id));
        if (old_par == nullptr) {
            std::string expl = "Parent OTT ID " + std::to_string(parent_id) + " is not an included taxon.";
            return bool_str_t{false, expl};
        }
        auto new_par = nd_ptr->get_parent();
        if (new_par != old_par) {
            old_par->remove_child(nd_ptr);
            new_par->add_child(nd_ptr);
        }
    }
    const auto & tr = get_new_tax_rec(oid, parent_id,
                                      name, rank, sourceinfo,
                                      uniqname, flags);
    if (name != nd_ptr->get_name() && name != nd_data.get_nonuniqname()) {
        remove_name_to_node_from_maps(nd_ptr->get_name(), nd_ptr);
        string cp{nd_data.get_nonuniqname()};
        remove_name_to_node_from_maps(cp, nd_ptr);
        auto mit = rt_data.name_to_node.find(nd_ptr->get_name());
        if (mit != rt_data.name_to_node.end() && mit->second == nd_ptr) {
            rt_data.name_to_node.erase(mit);
        }
        mit = rt_data.name_to_node.find(nd_data.get_nonuniqname());
        if (mit != rt_data.name_to_node.end() && mit->second == nd_ptr) {
            rt_data.name_to_node.erase(mit);
        }
    }
    reg_or_rereg_nd(nd_ptr, tr, tree);
    return bool_str_t{true, ""};
}

bool_str_t PatchableTaxonomy::add_new_taxon(OttId oid,
                                            OttId parent_id,
                                            const std::string & name,
                                            const std::string & rank,
                                            const std::string & sourceinfo,
                                            const std::string & uniqname,
                                            const std::string & flags,
                                            OttId * homonym_of) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
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
    auto itnit = included_taxon_from_id(oid);
    auto itrit = rt_data.id_to_record.find(oid);
    if (itnit != nullptr || itrit != rt_data.id_to_record.end()) {
        std::string expl = "OTT ID " + std::to_string(oid);
        expl += " is already used.";
        return bool_str_t{false, expl};
    }
    RTRichTaxNode * par_ptr = const_cast<RTRichTaxNode *>(included_taxon_from_id(parent_id));
    if (par_ptr == nullptr) {
        itrit = rt_data.id_to_record.find(parent_id);
        std::string expl = "Parent OTT ID " + std::to_string(parent_id);
        if (itrit != rt_data.id_to_record.end()) {
            expl += " refers to a filtered taxon.";
        } else {
            expl += " is unrecognized.";
        }
        return bool_str_t{false, expl};
    }
    if (uniqname.length() > 0) {
        return bool_str_t{false, "handling of uniqname not supported"};
    }

    const auto & tr = get_new_tax_rec(oid, parent_id,
                                      name, rank, sourceinfo,
                                      uniqname, flags);
    auto nnd = tree.create_child(par_ptr);
    reg_or_rereg_nd(nnd, tr, tree);
    return bool_str_t{true, ""};
}

void PatchableTaxonomy::reg_or_rereg_nd(RTRichTaxNode * nnd,
                                        const TaxonomyRecord & tr,
                                        RichTaxTree & tree) {
    auto nodeNamer = [](const auto&){return string();};
    populate_node_from_taxonomy_record(*nnd, tr, nodeNamer, tree);
    auto & rt_data = tree.get_data();
    rt_data.name_to_node[tr.name] = nnd;
    rt_data.id_to_node[nnd->get_ott_id()] =nnd;
}

TaxonomyRecord & PatchableTaxonomy::get_new_tax_rec(OttId oid,
                                                    OttId parent_id,
                                                    const std::string & name,
                                                    const std::string & rank,
                                                    const std::string & sourceinfo,
                                                    const std::string & uniqname,
                                                    const std::string & flags) {
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
    added_records.push_back(TaxonomyRecord(fake_line));
    return *(added_records.rbegin());
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
    for (auto tr_vsyn : rec_to_new_syn) {
        auto tr = tr_vsyn.first;
        for (auto new_syn : tr_vsyn.second) {
            sf << new_syn.name << sep 
                << tr->id << sep
                << sep // we don't retain the type on parsing
                << sep // we don't retain the uniqname on parsing
                << new_syn.source_string << sep << '\n';       
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
