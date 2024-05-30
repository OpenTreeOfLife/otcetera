#include <boost/algorithm/string/join.hpp>
#include <filesystem>
namespace fs = std::filesystem;
#include "otc/taxonomy/patching.h"
#include "otc/otc_base_includes.h"
#include "otc/ctrie/context_ctrie_db.h"
#include <boost/algorithm/string/join.hpp>

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
    for (auto node : iter_post_const(rich_tax_tree)) {
        assert(node != nullptr);
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

bool_str_t PatchableTaxonomy::delete_forward(OttId /* former_id */, OttId /*redirect_to_id*/) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();

    throw OTCError() << "delete_forward not implemented";
}

void possible_add_barren_flag_to_anc(RTRichTaxNode * nd_ptr) {
    if (nd_ptr == nullptr) {
        return;
    }
    int barren_bit = flag_from_string("barren");
    for (auto c : iter_child(*nd_ptr)) {
        RTRichTaxNodeData & c_data = c->get_data();
        if (!c_data.flags.test(barren_bit)) {
            return;
        }
    }
    RTRichTaxNodeData & nd_data = nd_ptr->get_data();
    nd_data.flags.set(barren_bit, true);
    possible_add_barren_flag_to_anc(nd_ptr->get_parent());
}

void possible_add_extinct_flag_to_anc(RTRichTaxNode * nd_ptr) {
    if (nd_ptr == nullptr) {
        return;
    }
    int ex_bit = flag_from_string("extinct");
    int ex_in_bit = flag_from_string("extinct_inherited");
    for (auto c : iter_child(*nd_ptr)) {
        RTRichTaxNodeData & c_data = c->get_data();
        if (! (c_data.flags.test(ex_bit) || c_data.flags.test(ex_in_bit))) {
            return;
        }
    }
    RTRichTaxNodeData & nd_data = nd_ptr->get_data();
    nd_data.flags.set(ex_in_bit, true);
    possible_add_extinct_flag_to_anc(nd_ptr->get_parent());
}

void detach_from_par_helper(RTRichTaxNode * nd_ptr, RTRichTaxNode * old_par) {
    old_par->remove_child(nd_ptr);
    if (!old_par->has_children()) {
        RTRichTaxNodeData & par_data = old_par->get_data();
        par_data.add_flags_from_string("barren");
        auto further_anc = old_par->get_parent();
        if (further_anc != nullptr) {
            possible_add_barren_flag_to_anc(further_anc);
        }
    } else {
        possible_add_barren_flag_to_anc(old_par);
    }
}

void flag_par_extinct_helper(RTRichTaxNode * par_nd) {
    if (!par_nd->has_children()) {
        throw OTCError() << "Expecting parent " << std::to_string(par_nd->get_ott_id()) << " to have children.";

    } else {
        possible_add_extinct_flag_to_anc(par_nd);
    }
}

bool_str_t PatchableTaxonomy::delete_taxon(OttId ott_id) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    auto nd_ptr = const_cast<RTRichTaxNode * >(included_taxon_from_id(ott_id));
    if (nd_ptr == nullptr) {
        std::string expl = "OTT ID " + std::to_string(ott_id) + " is not an included taxon.";
        return bool_str_t{false, expl};
    }
    RTRichTaxNodeData & nd_data = nd_ptr->get_data();
    auto old_par = nd_ptr->get_parent();
    if (nd_ptr->has_children()) {
        std::string expl = "OTT ID " + std::to_string(ott_id) + " has at least 1 child. Deletion of internal nodes is not currently supported.";
        return bool_str_t{false, expl};
    }
    rt_data.id_to_record.erase(ott_id);
    rt_data.id_to_node.erase(ott_id);
    detach_from_par_helper(nd_ptr, old_par);
    return bool_str_t{true, ""};
}

bool_str_t PatchableTaxonomy::delete_id_set(const OttIdSet & ott_id_set) {
    for (auto oid : ott_id_set) {
        auto sret = delete_taxon(oid);
        if (!sret.first) {
            return sret;
        }
    }
    return bool_str_t{true, ""};
}

bool_str_t PatchableTaxonomy::sink_taxon(OttId jr_oid, OttId sr_id) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    auto jr_nd_ptr = const_cast<RTRichTaxNode * >(included_taxon_from_id(jr_oid));
    if (jr_nd_ptr == nullptr) {
        std::string expl = "Junior OTT ID " + std::to_string(jr_oid) + " is not an included taxon.";
        return bool_str_t{false, expl};
    }
    auto sr_nd_ptr = const_cast<RTRichTaxNode * >(included_taxon_from_id(sr_id));
    if (sr_nd_ptr == nullptr) {
        std::string expl = "Senior OTT ID " + std::to_string(sr_id) + " is not an included taxon.";
        return bool_str_t{false, expl};
    }
    
    RTRichTaxNodeData & jr_nd_data = jr_nd_ptr->get_data();
    auto old_jr_par = jr_nd_ptr->get_parent();
    std::string jr_name = std::string{jr_nd_data.get_nonuniqname()};
    auto jr_src_vec = jr_nd_data.sourceinfoAsVec();
    std::set<std::string> jr_src_set{jr_src_vec.begin(), jr_src_vec.end()};

    auto jr_src_str = jr_nd_data.get_sources_as_fmt_str();
    auto dr = this->delete_taxon(jr_oid);
    if (!dr.first) {
        return dr;
    }
    auto asr = this->add_synonym(jr_name, sr_id, jr_src_str);
    if (!asr.first) {
        return asr;
    }
    RTRichTaxNodeData & sr_nd_data = sr_nd_ptr->get_data();
    auto sr_src_vec = sr_nd_data.sourceinfoAsVec();
    std::set<std::string> sr_src_set{sr_src_vec.begin(), sr_src_vec.end()};
    std::set<std::string> all_src = set_union_as_set(sr_src_set, jr_src_set);
    if (all_src != sr_src_set) {
        sr_nd_data.source_info = boost::algorithm::join(all_src, ",");
    }
    return this->add_forward(jr_oid, sr_id);
}

bool_str_t PatchableTaxonomy::append_prop_for_set(const OttIdSet & oids,
                                                  OttId parent,
                                                  const std::string & sourceinfo,
                                                  const tax_flags & flags) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    RTRichTaxNode * new_par_ptr = nullptr;
    if (parent != UINT_MAX && parent != 0) {
        new_par_ptr = const_cast<RTRichTaxNode * >(included_taxon_from_id(parent));
        if (new_par_ptr == nullptr) {
            std::string expl = "OTT ID for new parent" + std::to_string(parent) + " is not an included taxon.";
            return bool_str_t{false, expl};
        }
    }
    for (auto ott_id : oids) {
        auto nd_ptr = const_cast<RTRichTaxNode * >(included_taxon_from_id(ott_id));
        if (nd_ptr == nullptr) {
            std::string expl = "OTT ID " + std::to_string(ott_id) + " is not an included taxon.";
            return bool_str_t{false, expl};
        }
        if (new_par_ptr != nullptr) {
            auto old_par = nd_ptr->get_parent();
            if (new_par_ptr != old_par) {
                detach_from_par_helper(nd_ptr, old_par);
                new_par_ptr->add_child(nd_ptr);
            }
        }
        if (!sourceinfo.empty()) {
            throw OTCError() << "Bulk additions of sources not supported yet";
        }
        if (flags.any()) {
            RTRichTaxNodeData & nd_data = nd_ptr->get_data();
            nd_data.flags |=  flags;
            int extinct_bit = flag_from_string("extinct");
            if (flags.test(extinct_bit)) {
                auto old_par = nd_ptr->get_parent();
                flag_par_extinct_helper(old_par);
            }
        }
    }
    return bool_str_t{true, ""};
}

bool_str_t PatchableTaxonomy::edit_taxon(OttId oid,
                                         OttId parent_id,
                                         const std::string & name,
                                         const std::string & rank,
                                         const std::string & sourceinfo,
                                         const std::string & uniqname,
                                         const std::string & flags,
                                         bool flags_edited,
                                         OttId * /* homonym_of */) {
    auto & tree = this->get_mutable_tax_tree();
    auto & rt_data = tree.get_data();
    auto nd_ptr = const_cast<RTRichTaxNode * >(included_taxon_from_id(oid));
    if (nd_ptr == nullptr) {
        std::string expl = "OTT ID " + std::to_string(oid) + " is not an included taxon.";
        return bool_str_t{false, expl};
    }
    RTRichTaxNodeData & nd_data = nd_ptr->get_data();
    auto old_par = nd_ptr->get_parent();
    if (parent_id != UINT_MAX && parent_id != 0) {
        auto new_par = const_cast<RTRichTaxNode * >(included_taxon_from_id(parent_id));
        if (new_par == nullptr) {
            std::string expl = "Parent OTT ID " + std::to_string(parent_id) + " is not an included taxon.";
            return bool_str_t{false, expl};
        }
        if (new_par != old_par) {
            detach_from_par_helper(nd_ptr, old_par);
            new_par->add_child(nd_ptr);
        }
    } else if (old_par != nullptr) {
        parent_id = old_par->get_ott_id();
    }
    std::string fs = (flags_edited ? flags : flags_to_string(nd_data.get_flags()));
    std::string name_str = (name.empty() ? std::string{nd_data.get_nonuniqname()} : name);
    std::string rank_str = (rank.empty() ? nd_data.get_rank() : rank);
    std::string src_str = (sourceinfo.empty() ? nd_data.get_sources_as_fmt_str() : sourceinfo);
    std::string uname_str = uniqname;
    const auto & tr = get_new_tax_rec(oid, parent_id,
                                      name_str, rank_str, src_str,
                                      uname_str, fs);
    if (name_str != nd_ptr->get_name() && name_str != nd_data.get_nonuniqname()) {
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

    // 1. Check if the OTT ID is already used.
    auto itnit = included_taxon_from_id(oid);
    auto itrit = rt_data.id_to_record.find(oid);
    if (itnit != nullptr || itrit != rt_data.id_to_record.end()) {
        std::string expl = "OTT ID " + std::to_string(oid);
        expl += " is already used.";
        return bool_str_t{false, expl};
    }

    // 2. Check if the new taxon is a homonym of an existing taxon.
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

    // 3. Find the parent taxon by its OTT ID.
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

    // 4. Complain if the new taxon has a uniqname
    if (uniqname.length() > 0) {
        return bool_str_t{false, "handling of uniqname not supported"};
    }

    // 5. Create the new taxon record.
    const auto & tr = get_new_tax_rec(oid, parent_id,
                                      name, rank, sourceinfo,
                                      uniqname, flags);

    // 6. Create the new tree node
    auto nnd = tree.create_child(par_ptr);
    reg_or_rereg_nd(nnd, tr, tree);

    // 7. Update the fuzzy match databases
    if (auto f = get_fuzzy_matcher())
	f->add_key(name, oid, *this);

    return bool_str_t{true, ""};
}

void PatchableTaxonomy::reg_or_rereg_nd(RTRichTaxNode * nnd,
                                        const TaxonomyRecord & tr,
                                        RichTaxTree & tree) {
    auto nodeNamer = [](const auto&){return string();};
    populate_node_from_taxonomy_record(*nnd, tr, nodeNamer, tree);
    auto & rt_data = tree.get_data();
    rt_data.name_to_node[tr.name] = nnd;
    rt_data.id_to_node[nnd->get_ott_id()] = nnd;

    // Fill out depth field for node.
    nnd->get_data().depth = nnd->get_parent()->get_data().depth + 1;
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
