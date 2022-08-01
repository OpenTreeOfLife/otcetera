#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <bitset>
#include <regex>


#include "json.hpp"

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/diff_maker.h"
#include "otc/config_file.h"

INITIALIZE_EASYLOGGINGPP

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::unique_ptr;
using std::set;
using std::map;
using std::string_view;
using json = nlohmann::json;
using vec_strv_t = std::vector<std::string_view>;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;

variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;
    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("oldtaxonomy", value<string>(),"Filename for the old taxonomy")
        ("newtaxonomy", value<string>(),"Filename for the new taxonomy")
        ;

    options_description output("Output options");
    output.add_options()
        ("write-to-stdout","Primarily for debugging. Writes contents of taxonomy output to stdout. Only used if write-taxonomy is not used.")
        ;

    options_description visible;
    visible.add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("oldtaxonomy", 1);
    p.add("newtaxonomy", 1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-taxonomy-diff-maker <taxonomy-dir> [OPTIONS]\n"
                                                    "Read a taxonomy and edit JSON files",
                                                    visible, invisible, p);

    return vm;
}

using ndvec_t = std::vector<const RTRichTaxNode *>;
using id2name_t = std::unordered_map<OttId, string_view>;
using name2id_t = std::unordered_map<string_view, OttId>;
using nd2idset_t = std::unordered_map<const RTRichTaxNode *, OttIdSet>;
using idset2nd_vec_t = std::map<OttIdSet, ndvec_t >;

template<typename T>
bool fill_name_id_maps(const T & tree_data, id2name_t& id2name) {
    for (auto id_nd_pair : tree_data.id_to_node) {
        const auto & nd_ptr = id_nd_pair.second;
        const auto tax_id = id_nd_pair.first;
        if (nd_ptr == nullptr) {
            throw OTCError() << "Unexpected nullptr in id_to_node";
        }
        id2name[tax_id] = nd_ptr->get_name();
    }
    return true;
}

OttIdSet find_ids_with_same_names(const id2name_t & old_id2name, const id2name_t & new_id2name) {
    OttIdSet same_id_name;
    for (auto id_name_pair : old_id2name) {
        OttId tax_id = id_name_pair.first;
        auto new_it = new_id2name.find(tax_id);
        if (new_it != new_id2name.end()) {
            if (new_it->second == id_name_pair.second) {
                same_id_name.insert(tax_id);
            }
        }
    }
    return same_id_name;
}

nd2idset_t fill_term_des_id_set(const RichTaxTree & tree,
                                const OttIdSet & relevant_ids,
                                OttIdSet & terms_added) {
    nd2idset_t nd2idset;
    for (auto nd : iter_post_const(tree)) {
        const OttId tax_id = nd->get_ott_id();
        const bool is_relevant = contains(relevant_ids, tax_id);
        OttIdSet & dest = nd2idset[nd];
        // add a tip if it is not the ancestor of another relevant ID
        if (!nd->is_tip()) {
            for (auto child : iter_child_const(*nd)) {
                const auto & child_set = nd2idset.at(child);
                dest.insert(child_set.begin(), child_set.end());
            }
        }
        if (is_relevant && dest.empty()) {
            dest.insert(tax_id);
            terms_added.insert(tax_id);
        }
    }
    return nd2idset;
}

void partitionTaxonByTypeOfType(const RichTaxTree & tree,
                                OttIdSet & clade_ids,
                                OttIdSet & specimen_based) {
    std::vector<const RTRichTaxNode *> unk;
    for (auto nd : iter_pre_const(tree)) {
        auto & nd_data = nd->get_data();
        auto & rank = nd_data.rank;
        bool is_clade = false;
        const auto tax_id = nd->get_ott_id();
        assert(!contains(clade_ids, tax_id));
        assert(!contains(specimen_based, tax_id));
        if (rank == TaxonomicRank::RANK_NO_RANK && !nd->is_tip()) {
            unk.push_back(nd);
            continue;
        }
        if (rank < TaxonomicRank::RANK_SPECIES) {
            is_clade = true;
        }
        const auto & par = nd->get_parent();
        if (par) {
            if (is_clade) {
                if (contains(specimen_based, par->get_ott_id())) {
                    throw OTCError() << "Taxon with ID " << tax_id << " is higher taxon, but parent (" << par->get_ott_id() << ") is not.";
                }
            }
        }
        if (is_clade) {
            clade_ids.insert(tax_id);
        } else {
            specimen_based.insert(tax_id);
        }
    }
    while (unk.size() > 0) {
        std::vector<const RTRichTaxNode *> current;
        current.clear();
        std::swap(unk, current);
        const auto before_size = current.size();
        for (auto nd : current) {
            const auto tax_id = nd->get_ott_id();
            bool is_below_sp = false;
            for (auto anc : iter_anc_const(*nd)) {
                auto anc_id = anc->get_ott_id();
                if (contains(specimen_based, anc_id)) {
                    // LOG(DEBUG) << "Anc " << anc_id << " causing " << tax_id << " to be specimen_based.";
                    is_below_sp = true;
                    break;
                }
                if (contains(clade_ids, anc_id)) {
                    break;
                }
            }
            if (is_below_sp) {
                specimen_based.insert(tax_id);
                continue;
            }
            bool is_above_higher = false;
            for (auto c : iter_child_const(*nd)) {
                auto c_id = c->get_ott_id();
                if (contains(clade_ids, c_id)) {
                    // LOG(DEBUG) << "child " << c_id << " causing " << tax_id << " to be clade.";
                    is_above_higher = true;
                    break;
                }
            }
            if (is_above_higher) {
                clade_ids.insert(tax_id);
            } else {
                unk.push_back(nd);
            }
        }
        const auto after_size = unk.size();
        if (before_size == after_size) {
            // creeping from edges has stopped settling taxa
            //  call all of the rest as higher taxa
            for (auto nd : current) {
                clade_ids.insert(nd->get_ott_id());
            }
            unk.clear();
        }
        
    }
}

            

const RTRichTaxNode * find_specimen_based_root(const RTRichTaxNode *nd,
                                               const OttIdSet & spec_based_ids,
                                               OttIdSet & seen) {
    const RTRichTaxNode * par = nd->get_parent();
    const RTRichTaxNode * ret = nullptr;
    while (true) {
        if (!par) {
            ret = nd;
            break;
        }
        if (contains(spec_based_ids, par->get_ott_id())) {
            nd = par;
        } else {
            ret = nd;
            break;
        }
        par = nd->get_parent();
    }
    if (ret->is_tip()) {
        seen.insert(ret->get_ott_id());
    } else {
        for (auto nnd : iter_pre_n_const(ret)) {
            auto tax_id = nnd->get_ott_id();
            if (!contains(spec_based_ids, tax_id)) {
                throw OTCError() << "taxon " << nnd->get_name() << " (" << tax_id << ") not specimen_based, but in spec_based clade";
            }
            seen.insert(tax_id);
        }
    }
    return ret;
}

using id2nd_t = std::unordered_map<OttId, const RTRichTaxNode *>;
id2nd_t find_all_specimen_based_roots(const id2nd_t id2nd,
                                      const OttIdSet & specimen_based_ids,
                                      OttIdSet & seen) {
    id2nd_t sp_root;
    for (auto ott_id : specimen_based_ids) {
        if (contains(seen, ott_id)) {
            continue;
        }
        auto nd = id2nd.at(ott_id);
        auto root_ptr = find_specimen_based_root(nd, specimen_based_ids, seen);
        assert(root_ptr != nullptr);
        sp_root[root_ptr->get_ott_id()] = root_ptr;
    }
    return sp_root;
}


enum AlphaEditOp {
    NO_CHANGE = 0,
    CHANGED_ID = 1,
    CHANGED_NAME = 2,
    DELETED_SYN = 3,
    ADDED_SYN = 4,
    CHANGED_RANK = 5,
    DELETE_TAXON = 6,
    ADD_TAXON = 7
};

const vector<string> aeo2str = {"no change", "change id", "change name",
                                "delete synonym", "add synonym", 
                                "change_rank", "delete taxon", "add taxon"};


class AlphaEdit {
    public:
        AlphaEdit()
            :operation(AlphaEditOp::NO_CHANGE),
            first_rank(TaxonomicRank::RANK_NO_RANK),
            second_rank(TaxonomicRank::RANK_NO_RANK) {
        }

        AlphaEditOp operation;
        string first_str, second_str;
        OttId first_id, second_id;
        TaxonomicRank first_rank, second_rank;

        void add_to_json_array(json & jarr) const {
            if (operation == AlphaEditOp::NO_CHANGE) {
                return;
            } 
            json el;
            el["operation"] = aeo2str[operation];
            if (operation == AlphaEditOp::CHANGED_ID) {
                el["from"] = first_id;
                el["to"] = second_id;
            } else {
                el["taxon_id"] = first_id;
                if (operation == AlphaEditOp::CHANGED_NAME) {
                    el["from"] = first_str;
                    el["to"] = second_str;
                } else if (operation == AlphaEditOp::DELETED_SYN || operation == AlphaEditOp::ADDED_SYN) {
                    el["synonym"] = first_str;
                    if (!second_str.empty()) {
                        el["type"] = second_str;
                    }
                } else if (operation == AlphaEditOp::DELETE_TAXON) {
                    if (!first_str.empty()) {
                        el["name"] = first_str;
                    }
                } else if (operation == AlphaEditOp::CHANGED_RANK) {
                    el["from"] = rank_enum_to_name.at(first_rank);
                    el["to"] = rank_enum_to_name.at(second_rank);
                } else {
                    assert(operation == AlphaEditOp::ADD_TAXON);
                    el["name"] = first_str;
                    el["rank"] = rank_enum_to_name.at(first_rank);
                }
            }
            jarr.push_back(el);
        }
};

class TaxonomyDiffer {
    public:
    TaxonomyDiffer(const TaxonomyDiffMaker & old_tax, const TaxonomyDiffMaker & new_tax);
    void write(std::ostream & out ) const ;
    
    protected:
    void compare_specimen_based();
    void diagnose_root_spec_based_status(const RTRichTaxNode *old_spec_nd);
    void children_diagnose_root_spec_based_status(const RTRichTaxNode *old, const RTRichTaxNode * new_nd) {
    }
    void record_syn_diffs(const RTRichTaxNode *old_nd, const RTRichTaxNode *new_nd);

    AlphaEdit & new_alpha_edit();

    const RichTaxTree & old_tree;
    const RichTaxTree & new_tree;
    const RTRichTaxTreeData & old_td;
    const RTRichTaxTreeData & new_td;
    OttIdSet old_specimen_based_ids, new_specimen_based_ids;
    OttIdSet old_clade_ids, new_clade_ids;
    std::unordered_map<OttId, const RTRichTaxNode *> old_sp_root;
    std::unordered_map<string, ndvec_t> old_name_to_nds;
    std::unordered_map<string, ndvec_t> new_name_to_nds;

    std::vector<AlphaEdit> alphaTaxonomyEdits;
};

void TaxonomyDiffer::write(std::ostream & out ) const {
    json edits = json::array();
    for (const auto & el : alphaTaxonomyEdits) {
        el.add_to_json_array(edits);
    }
    json document;
    document["alpha"] = edits;
    out << document.dump(1) << std::endl;
}


TaxonomyDiffer::TaxonomyDiffer(const TaxonomyDiffMaker & old_tax,
                               const TaxonomyDiffMaker & new_tax)
  :old_tree(old_tax.get_tax_tree()),
   new_tree(new_tax.get_tax_tree()),
   old_td(old_tax.get_tax_tree().get_data()),
   new_td(new_tax.get_tax_tree().get_data()) {
    for (auto idnd_p : old_td.id_to_node) {
        const auto & nd_ptr = idnd_p.second;
        old_name_to_nds[nd_ptr->get_name()].push_back(nd_ptr);
    }
    for (auto idnd_p : new_td.id_to_node) {
        const auto & nd_ptr = idnd_p.second;
        new_name_to_nds[nd_ptr->get_name()].push_back(nd_ptr);
    }
    partitionTaxonByTypeOfType(old_tree, old_clade_ids, old_specimen_based_ids);
    LOG(DEBUG) << old_clade_ids.size() << " clades and " << old_specimen_based_ids.size() << " specimen_based ids.";
    partitionTaxonByTypeOfType(new_tree, new_clade_ids, new_specimen_based_ids);
    LOG(DEBUG) << new_clade_ids.size() << " clades and " << new_specimen_based_ids.size() << " specimen_based ids.";

    compare_specimen_based();
    
}

void TaxonomyDiffer::compare_specimen_based() {
    OttIdSet seen;
    old_sp_root = find_all_specimen_based_roots(old_td.id_to_node,
                                                old_specimen_based_ids, 
                                                seen);
    LOG(DEBUG) << "old tree had " << old_sp_root.size() << " roots of specimen_based taxa.";
    for (auto ott_id_root_p : old_sp_root) {
        diagnose_root_spec_based_status(ott_id_root_p.second);
    }
}

const AlphaEdit mtEdit;

inline AlphaEdit & TaxonomyDiffer::new_alpha_edit() {
    alphaTaxonomyEdits.push_back(mtEdit);
    return *(alphaTaxonomyEdits.rbegin());
}

void TaxonomyDiffer::record_syn_diffs(const RTRichTaxNode *old_nd, const RTRichTaxNode *new_nd) {
    if (new_nd == nullptr) {
        assert(old_nd != nullptr);
        const RTRichTaxNodeData & old_nd_data = old_nd->get_data();
        for (const auto js : old_nd_data.junior_synonyms) {
            auto & edit = new_alpha_edit();
            edit.operation = AlphaEditOp::DELETED_SYN;
            edit.first_id = old_nd->get_ott_id();
            edit.first_str = js->name;
            if (!js->source_string.empty()) {
                edit.second_str = js->source_string;
            }
        }
        return;
    }
    if (old_nd == nullptr) {
        assert(new_nd != nullptr);
        const RTRichTaxNodeData & new_nd_data = new_nd->get_data();
        for (const auto js : new_nd_data.junior_synonyms) {
            auto & edit = new_alpha_edit();
            edit.operation = AlphaEditOp::ADDED_SYN;
            edit.first_id = new_nd->get_ott_id();
            edit.first_str = js->name;
            if (!js->source_string.empty()) {
                edit.second_str = js->source_string;
            }
        }
        return;
    }
    map<string, string> oldpairs;
    map<string, string> newpairs;
    const RTRichTaxNodeData & old_nd_data = old_nd->get_data();
    for (const auto js : old_nd_data.junior_synonyms) {
        oldpairs[js->name] = js->source_string;
    }
    const RTRichTaxNodeData & new_nd_data = new_nd->get_data();
    for (const auto js : new_nd_data.junior_synonyms) {
        newpairs[js->name] = js->source_string;
    }
    for (auto op : oldpairs) {
        auto nIt = newpairs.find(op.first);
        if (nIt == newpairs.end()) {
            auto & edit = new_alpha_edit();
            edit.operation = AlphaEditOp::DELETED_SYN;
            edit.first_id = new_nd->get_ott_id();
            edit.first_str = op.first;
            edit.second_str = op.second;
            continue;
        }
        if (nIt->second != op.second) {
            auto & edit = new_alpha_edit();
            edit.operation = AlphaEditOp::DELETED_SYN;
            edit.first_id = new_nd->get_ott_id();
            edit.first_str = op.first;
            edit.second_str = op.second;
            edit = new_alpha_edit();
            edit.operation = AlphaEditOp::ADDED_SYN;
            edit.first_id = new_nd->get_ott_id();
            edit.first_str = nIt->first;
            edit.second_str = nIt->second;
        }
    }
    for (auto np : newpairs) {
        auto oIt = oldpairs.find(np.first);
        if (oIt == oldpairs.end()) {
            auto & edit = new_alpha_edit();
            edit.operation = AlphaEditOp::ADDED_SYN;
            edit.first_id = new_nd->get_ott_id();
            edit.first_str = np.first;
            edit.second_str = np.second;
        }
    }
}

void TaxonomyDiffer::diagnose_root_spec_based_status(const RTRichTaxNode *old_spec_nd) {
    const auto & new_i2nd = new_td.id_to_node;
    const OttId old_root_id = old_spec_nd->get_ott_id();
    const auto new_el = new_i2nd.find(old_root_id);
    const RTRichTaxNode * new_cmp_nd = nullptr;

    if (new_el == new_i2nd.end()) {
        const auto old_name = old_spec_nd->get_name();
        const auto new_by_name = new_name_to_nds.find(old_name);
        if (new_by_name == new_name_to_nds.end()) {
            children_diagnose_root_spec_based_status(old_spec_nd, nullptr);
            // need to keep this before the DELETE_TAXON here
            record_syn_diffs(old_spec_nd, new_cmp_nd);
            auto & edit = new_alpha_edit();
            edit.operation = AlphaEditOp::DELETE_TAXON;
            edit.first_id = old_root_id;
            edit.first_str = old_name;
        } else {
            const ndvec_t newnodes = new_by_name->second;
            if (newnodes.size() != 1) {
                throw OTCError() << "Old specimen_based root " << old_name << " (" << old_root_id << ") not found by id. and maps to multiple names";
            }
            new_cmp_nd = newnodes[0];
            children_diagnose_root_spec_based_status(old_spec_nd, new_cmp_nd);
            record_syn_diffs(old_spec_nd, new_cmp_nd); 
            auto & edit = new_alpha_edit();
            edit.operation = AlphaEditOp::CHANGED_NAME;
            edit.first_str = old_name;
            edit.second_str = new_cmp_nd->get_name();
        }
    } else {
        new_cmp_nd = new_el->second;
        children_diagnose_root_spec_based_status(old_spec_nd, new_el->second);
        record_syn_diffs(old_spec_nd, new_cmp_nd);
    }
}

bool diff_from_taxonomies(std::ostream & out,
                          const TaxonomyDiffMaker & old_tax,
                          const TaxonomyDiffMaker & new_tax) {
    TaxonomyDiffer tax_dif(old_tax, new_tax);
    tax_dif.write(out);
    // // I. we look at the subset of taxa that have the same ID and name between versions
    // //   I.A - find the set of IDs with the same name
    // id2name_t old_id2name, new_id2name;
    // LOG(DEBUG) << "old_td.name_to_node.size() = " << old_td.name_to_node.size() ;
    // LOG(DEBUG) << "new_td.name_to_node.size() = " << new_td.name_to_node.size() ;
    // fill_name_id_maps(old_td, old_id2name);
    // fill_name_id_maps(new_td, new_id2name);
    // LOG(DEBUG) << old_id2name.size() << " " << new_id2name.size() << std::endl;
    // const OttIdSet same_id_name = find_ids_with_same_names(old_id2name, new_id2name);
    // LOG(DEBUG) << same_id_name.size() << " IDs with the same name between versions";
    // //   I.B - further restrict this (as "culled") to the set of those IDs that are
    // //      terminal when only the same_id_name IDs are relevant.
    // OttIdSet ota, nta;
    // fill_term_des_id_set(old_tree, same_id_name, ota);
    // fill_term_des_id_set(new_tree, same_id_name, nta);
    // const auto & tmp = intersection_of_sets(ota, nta);
    // OttIdSet culled = intersection_of_sets(same_id_name, tmp);
    // ota.clear();
    // nta.clear();
    // //   I.C - get the mappings from ID to relevant set of IDs for each tree.
    // const nd2idset_t old_nd2ids = fill_term_des_id_set(old_tree, culled, ota);
    // const nd2idset_t new_nd2ids = fill_term_des_id_set(new_tree, culled, nta);
    // assert(ota == nta);
    // assert(culled == nta);
    // LOG(DEBUG) << culled.size() << " = culled.size()";
    // //   I.D - look for phylorefs that are new or deleted wrt to the culled set of IDs
    // idset2nd_vec_t old_idset2ndvec, new_idset2ndvec;
    // for (auto nidsp: old_nd2ids) {
    //     old_idset2ndvec[nidsp.second].push_back(nidsp.first);
    // }
    // LOG(DEBUG) << old_idset2ndvec.size() << " = old_idset2ndvec.size()";
    // for (auto nidsp: new_nd2ids) {
    //     new_idset2ndvec[nidsp.second].push_back(nidsp.first);
    // }
    // LOG(DEBUG) << new_idset2ndvec.size() << " = new_idset2ndvec.size()";
    // for (auto idsnp: new_idset2ndvec) {
    //     auto & idset = idsnp.first;
    //     if (!contains(old_idset2ndvec, idset)) {
    //         out << "INSERT MRCA(";
    //         write_tax_id_set(out, "", idset, ",");
    //         out << ")" << std::endl;
    //     }
    // }
    
    return true;
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc, argv);
        std::ostream & out = std::cout;
        if (!args.count("oldtaxonomy")) {
            cerr << "oldtaxonomy expected as first unnamed argument\n";
            return 1;
        }
        if (!args.count("newtaxonomy")) {
            cerr << "newtaxonomy expected as second unnamed argument\n";
            return 1;
        }
        string otd = args["oldtaxonomy"].as<string>();
        string ntd = args["newtaxonomy"].as<string>();
        OttId keep_root = -1;
        bitset<32> cleaning_flags = 0;
        LOG(INFO) << "loading old taxonomy\n";
        Taxonomy::tolerate_synonyms_to_unknown_id = true;
        TaxonomyDiffMaker otaxonomy = {otd, cleaning_flags, keep_root};
        LOG(INFO) << "loading new taxonomy\n";
        TaxonomyDiffMaker ntaxonomy = {ntd, cleaning_flags, keep_root};
        diff_from_taxonomies(out, otaxonomy, ntaxonomy);
        
    } catch (std::exception& e) {
        cerr << "otc-taxonomy-diff-maker: Error! " << e.what() << std::endl;
        return 1;
    }
}

// 1. Write a parser to read the lines faster
// 2. Avoid memory allocation -- by mmapping the taxonomy file?
// 3. Convert the flags into a bitmask
// 4. Should the Rank be a converted to an integer?
// 5. Can we assign OTT IDs to internal nodes of a tree while accounting for Incertae Sedis taxa?
// * What are the triplet-inference rules for the Incertae Sedis problem?

// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

