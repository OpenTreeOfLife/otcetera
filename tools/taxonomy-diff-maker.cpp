#include "otc/taxonomy/taxonomy-diff.h"
#include <boost/program_options.hpp>
#include "otc/taxonomy/diff_maker.h"

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
using id2nd_t = std::unordered_map<OttId, const RTRichTaxNode *>;


using id2grouping_t = map<OttId, Grouping>;
using name2grouping_t = map<string, vector<Grouping *> >;

class TaxonomyDiffer {
    public:
    TaxonomyDiffer(TaxonomyDiffMaker & old_tax, TaxonomyDiffMaker & new_tax);
    void write(std::ostream & out ) const ;
    
    protected:
    void compare_specimen_based();
    void compare_higher_taxa();
    //void diagnose_fate_of_groupings(const RTRichTaxNode *old_spec_nd);
    void diagnose_old_spec_based_fate(const RTRichTaxNode *old_spec_nd, bool top_level);
    void find_pair_for_new(const RTRichTaxNode *new_nd,
                           id2grouping_t & old_by_id,
                           name2grouping_t & old_by_name,
                           id2grouping_t & new_by_id);
    //void diagnose_new_spec_based_status(const RTRichTaxNode *old_spec_nd, bool top_level);
    void children_diagnose_old_spec_based_fate(const RTRichTaxNode *old,
                                               const RTRichTaxNode * new_nd,
                                               bool top_level);
    void record_syn_diffs(const RTRichTaxNode *old_nd, const RTRichTaxNode *new_nd);
    void record_tax_edits_for_match(const RTRichTaxNode *old_nd, const RTRichTaxNode * new_nd);

    AlphaEdit & new_alpha_edit();
    AlphaGroupEdit & new_alpha_group_edit();
    AlphaGroupEdit & new_higher_edit();

    RichTaxTree & old_tree;
    RichTaxTree & new_tree;
    const RTRichTaxTreeData & old_td;
    const RTRichTaxTreeData & new_td;
    OttIdSet old_specimen_based_ids, new_specimen_based_ids;
    OttIdSet old_clade_ids, new_clade_ids;
    std::unordered_map<OttId, const RTRichTaxNode *> old_sp_root;
    std::unordered_map<OttId, const RTRichTaxNode *> new_sp_root;
    std::unordered_map<string, ndvec_t> old_name_to_nds;
    std::unordered_map<string, ndvec_t> new_name_to_nds;

    // These 3 fields are filled starting with children_diagnose_old_spec_based_fate
    OttIdSet retained_spec_ids;  // found in both - same ID
    map<OttId, OttId> mapped_spec_ids;  // found in both. old ID -> new ID
    map<OttId, OttId> revmapped_spec_ids; //              new ID -> old ID
    OttIdSet deleted_spec_ids; // only in old
    OttIdSet new_spec_ids;     // only in new.

    // used in diagnose_fate_of_groupings
    OttIdSet new_handled, old_handled;
    
    std::vector<AlphaEdit> alphaTaxonomyEdits;
    std::vector<AlphaGroupEdit> alphaGroupEdits;
    std::vector<AlphaGroupEdit> higherGroupEdits;
};

const AlphaEdit mtEdit;
const AlphaGroupEdit mtAGEdit;

inline AlphaEdit & TaxonomyDiffer::new_alpha_edit() {
    alphaTaxonomyEdits.push_back(mtEdit);
    return *(alphaTaxonomyEdits.rbegin());
}

inline AlphaGroupEdit & TaxonomyDiffer::new_alpha_group_edit() {
    alphaGroupEdits.push_back(mtAGEdit);
    return *(alphaGroupEdits.rbegin());
}

inline AlphaGroupEdit & TaxonomyDiffer::new_higher_edit() {
    higherGroupEdits.push_back(mtAGEdit);
    return *(higherGroupEdits.rbegin());
}



void TaxonomyDiffer::write(std::ostream & out ) const {
    json edits = json::array();
    for (const auto & el : alphaTaxonomyEdits) {
        el.add_to_json_array(edits);
    }
    json gredits = json::array();
    for (const auto & grel : alphaGroupEdits) {
        grel.add_to_json_array(gredits);
    }
    json hedits = json::array();
    for (const auto & hrel : higherGroupEdits) {
        hrel.add_to_json_array(hedits);
    }
    json document;
    document["alpha"] = edits;
    document["alpha_groups"] = gredits;
    document["higher_taxa"] = hedits;
    out << document.dump(0) << std::endl;
}


TaxonomyDiffer::TaxonomyDiffer(TaxonomyDiffMaker & old_tax,
                               TaxonomyDiffMaker & new_tax)
  :old_tree(const_cast<RichTaxTree &>(old_tax.get_tax_tree())),
   new_tree(const_cast<RichTaxTree &>(new_tax.get_tax_tree())),
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
    compare_higher_taxa();
    
}

void fill_from_children(const RTRichTaxNode * inner_nd,
                        OttIdSet & uniq,
                        OttIdSet & retained,
                        const nd2idset_t & nd2uniq_ids,
                        const nd2idset_t & nd2des_ids) {
    for (auto c : iter_child_const(*inner_nd)) {
        const auto cuniq = nd2uniq_ids.find(c);
        if (cuniq != nd2uniq_ids.end()) {
            uniq.insert(cuniq->second.begin(), cuniq->second.end());
        }
        const auto cretained = nd2des_ids.find(c);
        if (cretained != nd2des_ids.end()) {
            retained.insert(cretained->second.begin(), cretained->second.end());
        }
    }
}

void TaxonomyDiffer::compare_higher_taxa() {
    map<std::uint32_t, const RTRichTaxNode *> old_by_trav, new_by_trav;
    set_traversal_entry_exit(old_tree);
    set_traversal_entry_exit(new_tree);
    nd2idset_t nd2add_ids;
    nd2idset_t nd2des_ids;
    for (auto nsrIt : new_sp_root) {
        const auto sp_id = nsrIt.first;
        auto sp_nd = nsrIt.second;
        if (contains(new_spec_ids, sp_id)) {
            nd2add_ids[sp_nd].insert(sp_id);
        } else {
            nd2des_ids[sp_nd].insert(sp_id);
        }
        for (auto anc : iter_anc_const(*sp_nd)) {
            const auto trav_enter = anc->get_data().trav_enter;
            if (contains(new_by_trav, trav_enter)) {
                break;
            }
            new_by_trav[trav_enter] = anc;
        }
    }
    LOG(DEBUG) << "new_by_trav.size() = " << new_by_trav.size() << " from " << new_by_trav.begin()->first << " " << new_by_trav.begin()->second->get_name() << " to " << new_by_trav.rbegin()->first << " " << new_by_trav.rbegin()->second->get_name();
    nd2idset_t nd2del_ids;
    nd2idset_t nd2des_ids_in_new;
    // fill in del and des fields for the "tips" of the tree (the species, presumably)
    for (auto osrIt : old_sp_root) {
        const auto sp_id = osrIt.first;
        auto sp_nd = osrIt.second;
        if (contains(deleted_spec_ids, sp_id)) {
            nd2del_ids[sp_nd].insert(sp_id);
        } else if (contains(retained_spec_ids, sp_id)) {
            nd2des_ids_in_new[sp_nd].insert(sp_id);
        } else {
            nd2des_ids_in_new[sp_nd].insert(mapped_spec_ids.at(sp_id));
        }
        for (auto anc : iter_anc_const(*sp_nd)) {
            const auto trav_enter = anc->get_data().trav_enter;
            if (contains(old_by_trav, trav_enter)) {
                break;
            }
            old_by_trav[trav_enter] = anc;
        }
    }
    LOG(DEBUG) << "old_by_trav.size() = " << old_by_trav.size() << " from " << old_by_trav.begin()->first << " " << old_by_trav.begin()->second->get_name() << " to " << old_by_trav.rbegin()->first << " " << old_by_trav.rbegin()->second->get_name();
    const auto & old_i2nd = old_td.id_to_node;
    const auto & new_i2nd = new_td.id_to_node;

    for (auto byIt = old_by_trav.rbegin(); byIt != old_by_trav.rend(); ++byIt) {
        auto inner_nd = byIt->second;
        OttIdSet & del = nd2del_ids[inner_nd];
        OttIdSet & in_new = nd2des_ids_in_new[inner_nd];
        auto tax_id = inner_nd->get_ott_id();
        auto ni2nIt = new_i2nd.find(tax_id);
        if (ni2nIt == new_i2nd.end()) {
            fill_from_children(inner_nd, del, in_new, nd2del_ids, nd2des_ids_in_new);
            auto & edit = new_higher_edit();
            edit.operation = AlphaGroupEditOp::DELETED_GROUPING;
            edit.first_id = tax_id;
            edit.first_str = inner_nd->get_name();
        } else {
            in_new.insert(tax_id);
        }
    }

    for (auto byIt = new_by_trav.rbegin(); byIt != new_by_trav.rend(); ++byIt) {
        auto inner_nd = byIt->second;
        OttIdSet & added = nd2add_ids[inner_nd];
        OttIdSet & retained = nd2des_ids[inner_nd];
        auto tax_id = inner_nd->get_ott_id();
        auto oi2nIt = old_i2nd.find(tax_id);
        if (oi2nIt == old_i2nd.end()) {
            fill_from_children(inner_nd, added, retained, nd2add_ids, nd2des_ids);
            auto & edit = new_higher_edit();
            edit.operation = AlphaGroupEditOp::NEW_GROUPING;
            edit.first_id = tax_id;
            edit.first_str = inner_nd->get_name();
            edit.addedIds = retained;
            edit.addedIds.insert(added.begin(), added.end());
        } else {
            retained.insert(tax_id);
            OttIdSet new_add, new_retained;
            fill_from_children(inner_nd, new_add, new_retained, nd2add_ids, nd2des_ids);
            auto old_nd = oi2nIt->second;
            OttIdSet old_del, old_retained;
            fill_from_children(old_nd, old_del, old_retained, nd2del_ids, nd2des_ids_in_new);
            if (inner_nd->get_name() != old_nd->get_name()) {
                auto & edit = new_higher_edit();
                edit.operation = AlphaGroupEditOp::GR_CHANGED_NAME;
                edit.first_id = tax_id;
                edit.first_str = old_nd->get_name();
                edit.second_str = inner_nd->get_name();
            }
            OttIdSet ret_but_del = set_difference_as_set(old_retained, new_retained);
            OttIdSet ret_but_add = set_difference_as_set(new_retained, old_retained);
            old_del.insert(ret_but_del.begin(), ret_but_del.end());
            new_add.insert(ret_but_add.begin(), ret_but_add.end());
            if (!new_add.empty()) {
                auto & edit = new_higher_edit();
                edit.first_id = tax_id;
                if (!old_del.empty()) {
                    edit.operation = AlphaGroupEditOp::ADD_DEL_TAXA;
                    edit.delIds = old_del;
                } else {
                    edit.operation = AlphaGroupEditOp::ADD_TAXA;
                }
                edit.addedIds = new_add;
            } else if (!old_del.empty()) {
                auto & edit = new_higher_edit();
                edit.first_id = tax_id;
                edit.operation = AlphaGroupEditOp::DEL_TAXA;
                edit.delIds = old_del;
            }
        }
    }
}


void TaxonomyDiffer::compare_specimen_based() {
    OttIdSet seen;
    const auto & old_i2nd = old_td.id_to_node;
    old_sp_root = find_all_specimen_based_roots(old_i2nd, old_specimen_based_ids, seen);
    OttIdSet nseen;
    const auto & new_i2nd = new_td.id_to_node;
    LOG(DEBUG) << "old tree had " << old_sp_root.size() << " roots of specimen_based taxa.";
    // this loop and recursive call will record all deleted, remapped and retained.
    for (auto ott_id_root_p : old_sp_root) {
        diagnose_old_spec_based_fate(ott_id_root_p.second, true);
    }
    LOG(DEBUG) << "retained_spec_ids.size() = " << retained_spec_ids.size();
    LOG(DEBUG) << "mapped_spec_ids.size() = " << mapped_spec_ids.size();
    LOG(DEBUG) << "revmapped_spec_ids.size() = " << revmapped_spec_ids.size();
    LOG(DEBUG) << "deleted_spec_ids.size() = " << deleted_spec_ids.size();
    // Now we can figure out the added specimen-based taxa
    for (auto nid : new_specimen_based_ids) {
        auto nd = new_i2nd.at(nid);
        assert(nd != nullptr);
        auto new_tax_id = nd->get_ott_id();
        if (contains(retained_spec_ids, new_tax_id) || 
            contains(revmapped_spec_ids, new_tax_id)) {
            continue;
        } 
        auto & edit = new_alpha_edit();
        edit.operation = AlphaEditOp::ADD_TAXON;
        edit.first_id = nid;
        edit.first_str = nd->get_name();
        auto & nd_data = nd->get_data();
        edit.first_rank = nd_data.rank;
        edit.first_flags = nd_data.flags;
        record_syn_diffs(nullptr, nd);
        new_spec_ids.insert(new_tax_id);
    }
    LOG(DEBUG) << "new_spec_ids.size() = " << new_spec_ids.size();
    id2grouping_t old_groups;
    name2grouping_t groups_by_name;
    for (auto ott_id_root_p : old_sp_root) {
        auto old_nd = ott_id_root_p.second;
        if (old_nd->is_tip()) {
            continue;
        }
        map<const RTRichTaxNode *, OttIdSet> ndToSet;
        for (auto desnd : iter_post_n_const(*old_nd)) {
            if (old_nd != desnd) {
                const RTRichTaxNode * par = desnd->get_parent();
                if (desnd->is_tip()) {
                    ndToSet[par].insert(desnd->get_ott_id());
                    continue;
                }
                const OttIdSet & self_set = ndToSet.at(desnd);
                ndToSet[par].insert(self_set.begin(), self_set.end());
            }
            Grouping & group = old_groups[desnd->get_ott_id()];
            const OttIdSet & curr_set = ndToSet.at(desnd);
            group.name = desnd->get_name();
            groups_by_name[group.name].push_back(&group); // register in map by name
            group.tax_id = desnd->get_ott_id();
            for (auto oid : curr_set) {
                if (contains(deleted_spec_ids, oid)) {
                    //NO-OP 
                } else if (contains(retained_spec_ids, oid)) {
                    group.shared_ids.insert(oid);
                } else {
                    try {
                        group.shared_ids.insert(mapped_spec_ids.at(oid));    
                    } catch (...) {
                        throw OTCError() << "Could not find ID " << oid << " (a des. of " << group.tax_id << ") in spec_id containers.";
                    }
                }
            }
        }
    }
    LOG(DEBUG) << "old_groups.size() = " << old_groups.size();
    LOG(DEBUG) << "groups_by_name.size() = " << groups_by_name.size();

    id2grouping_t tough_new_groups;
    new_sp_root = find_all_specimen_based_roots(new_i2nd, new_specimen_based_ids, nseen);
    for (auto ott_id_root_p : new_sp_root) {
        if (!ott_id_root_p.second->is_tip()) {
            find_pair_for_new(ott_id_root_p.second, old_groups, groups_by_name, tough_new_groups);
        }
    }
    id2grouping_t unpaired_old_groups;
    for (auto ogIt : old_groups) {
        if (!ogIt.second.paired) {
            unpaired_old_groups[ogIt.first] = ogIt.second;
        }
    }
    old_groups.clear();
    groups_by_name.clear();
    LOG(DEBUG) << "unpaired_old_groups.size() = " << unpaired_old_groups.size();
    LOG(DEBUG) << "tough_new_groups.size() = " << tough_new_groups.size();
    for (auto uogIt : unpaired_old_groups) {
        auto & agedit = new_alpha_group_edit();
        agedit.operation = AlphaGroupEditOp::DELETED_GROUPING;
        agedit.first_id = uogIt.first;
    }
    for (auto tngIt : tough_new_groups) {
        auto & agedit = new_alpha_group_edit();
        agedit.operation = AlphaGroupEditOp::NEW_GROUPING;
        agedit.first_id = tngIt.first;
        Grouping & grouping = tngIt.second;
        OttIdSet added_ids = grouping.shared_ids;
        added_ids.insert(grouping.new_ids.begin(), grouping.new_ids.end());
        agedit.addedIds = added_ids;
    }
    // now we characterize the structure of groupings below the species level.
    // using pair_nd_t = std::pair<const RTRichTaxNode *, const RTRichTaxNode *>;
    // vector<pair_nd_t> matched;
    // matched.reserve(retained_spec_ids.size() + mapped_spec_ids.size());
    // vector<const RTRichTaxNode *old_nd> old_only_intern, new_only_intern;
    // for (auto rid : retained_spec_ids) {
    //     auto old_nd = old_i2nd.at(rid);
    //     auto new_nd = new_i2nd.at(rid);
    //     if (old_nd->is_tip()) {
    //         if (!new_nd->is_tip()) {
    //             new_only_intern.push_back(new_nd);
    //         }
            
    //     } else {
    //         if (new_nd->is_tip()) {
    //             old_only_intern.push_back(old_only_intern);
    //         } else {
    //             matched.push_back(pair_nd_t{old_nd, new_nd});
    //         }
    //     }
        
    // }
    // for (auto mIt : mapped_spec_ids) {
    //     auto old_nd = old_i2nd.at(mIt.first);
    //     auto new_nd = new_i2nd.at(mIt.second);
    //     matched.push_back(pair_nd_t{old_nd, new_nd});
    // }

}


void TaxonomyDiffer::find_pair_for_new(const RTRichTaxNode *new_nd,
                                       id2grouping_t & old_gr_by_id,
                                       name2grouping_t & old_gr_by_name,
                                       id2grouping_t & nonobvious
                                       ) {
    map<const RTRichTaxNode *, OttIdSet> ndToSet;
    id2grouping_t new_groups;
    for (auto desnd : iter_post_n_const(*new_nd)) {
        if (new_nd != desnd) {
            const RTRichTaxNode * par = desnd->get_parent();
            if (desnd->is_tip()) {
                ndToSet[par].insert(desnd->get_ott_id());
                continue;
            }
            const OttIdSet & self_set = ndToSet.at(desnd);
            ndToSet[par].insert(self_set.begin(), self_set.end());
        }
        Grouping & group = new_groups[desnd->get_ott_id()];
        const OttIdSet & curr_set = ndToSet.at(desnd);
        group.name = desnd->get_name();
        group.tax_id = desnd->get_ott_id();
        for (auto oid : curr_set) {
            if (contains(new_spec_ids, oid)) {
                group.new_ids.insert(oid);
            } else  {
                assert(contains(retained_spec_ids, oid) || contains(revmapped_spec_ids, oid));
                group.shared_ids.insert(oid);
            }
        }
    }

    for (auto grIt : new_groups) {
        OttId new_id = grIt.first;
        Grouping & grouping = grIt.second;
        auto rmIt = revmapped_spec_ids.find(new_id);
        const OttId correspond_id = (rmIt == revmapped_spec_ids.end() ? new_id : rmIt->second);
        auto old_grIt = old_gr_by_id.find(correspond_id);
        Grouping * old_grouping = nullptr;
        if (old_grIt == old_gr_by_id.end() || old_grIt->second.paired) {
            auto old_bynIt = old_gr_by_name.find(grouping.name);
            if (old_bynIt != old_gr_by_name.end()) {
                throw OTCError() << "not implemented - grouping found by name, but not ID";
            }
        } else {
            old_grouping = &(old_grIt->second);
        }
        if (old_grouping == nullptr) {
            nonobvious[new_id] = grouping;
        } else {
            grouping.paired = true;
            old_grouping->paired = true;
            if (correspond_id != new_id) {
                auto & agedit = new_alpha_group_edit();
                agedit.operation = AlphaGroupEditOp::GR_CHANGED_ID;
                agedit.first_id = correspond_id;
                agedit.second_id = new_id;
            }
            // if (grouping.name != old_grouping->name) {
            //     auto & agedit = new_alpha_group_edit();
            //     agedit.first_id = new_id;
            //     agedit.operation = AlphaGroupEditOp::GR_CHANGED_NAME;
            //     agedit.first_str = old_grouping->name;
            //     agedit.second_str = grouping.name;
            // }
            OttIdSet del_ids = set_difference_as_set(old_grouping->shared_ids, grouping.shared_ids);
            OttIdSet added_ids = set_difference_as_set(grouping.shared_ids, old_grouping->shared_ids);
            added_ids.insert(grouping.new_ids.begin(), grouping.new_ids.end());
            if (!added_ids.empty()) {
                auto & agedit = new_alpha_group_edit();
                agedit.first_id = new_id;
                if (!del_ids.empty()) {
                    agedit.operation = AlphaGroupEditOp::ADD_DEL_TAXA;
                    agedit.delIds = del_ids;
                } else {
                    agedit.operation = AlphaGroupEditOp::ADD_TAXA;
                }
                agedit.addedIds = added_ids;
            } else if (!del_ids.empty()) {
                auto & agedit = new_alpha_group_edit();
                agedit.first_id = new_id;
                agedit.operation = AlphaGroupEditOp::DEL_TAXA;
                agedit.delIds = del_ids;
            }
        }
    }

    // for (auto oid : retained_spec_ids) {
    //     find_pair_for_new(oid, &oid, old_groups, groups_by_name);
    // }
    // for (auto rmIt : revmapped_spec_ids) {
    //     OttId old_id = rmIt.second;
    //     find_pair_for_new(rmIt.first, &old_id, old_groups, groups_by_name);
    // }
    // for (auto oid : new_spec_ids) {
    //     find_pair_for_new(oid, nullptr, old_groups, groups_by_name);
    // }
    
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

void TaxonomyDiffer::children_diagnose_old_spec_based_fate(const RTRichTaxNode * old_nd,
                                                              const RTRichTaxNode * new_nd,
                                                              bool top_level) {
    assert(old_nd != nullptr);
    const auto old_tax_id = old_nd->get_ott_id();
    if (new_nd == nullptr) {
        deleted_spec_ids.insert(old_tax_id);
    } else {
        const auto new_tax_id = new_nd->get_ott_id();
        if (new_nd->get_ott_id() == old_tax_id) {
            retained_spec_ids.insert(old_tax_id);
        } else {
            mapped_spec_ids[old_tax_id] = new_tax_id;
            revmapped_spec_ids[new_tax_id] = old_tax_id;
        }
    }

    // indirect recursion
    for (auto child : iter_child_const(*old_nd)) {
        diagnose_old_spec_based_fate(child, false);
    }
    if (!top_level) {
        return;
    }
    // for (auto child : iter_child_const(*new_nd)) {
    //     diagnose_new_spec_based_status(child);
    // }
}

// void TaxonomyDiffer::diagnose_new_spec_based_status(const RTRichTaxNode *new_spec_nd) {
// }

void TaxonomyDiffer::diagnose_old_spec_based_fate(const RTRichTaxNode *old_spec_nd,
                                                  bool top_level) {
    const auto & new_i2nd = new_td.id_to_node;
    const OttId old_root_id = old_spec_nd->get_ott_id();
    const auto new_el = new_i2nd.find(old_root_id);
    const RTRichTaxNode * new_cmp_nd = nullptr;

    if (new_el == new_i2nd.end()) {
        const auto old_name = old_spec_nd->get_name();
        const auto new_by_name = new_name_to_nds.find(old_name);
        if (new_by_name == new_name_to_nds.end()) {
            children_diagnose_old_spec_based_fate(old_spec_nd, nullptr, top_level);
            // don't need to worry about deleting synonyms of a deleted  taxon
            //record_syn_diffs(old_spec_nd, new_cmp_nd);
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
            children_diagnose_old_spec_based_fate(old_spec_nd, new_cmp_nd, top_level);
            record_tax_edits_for_match(old_spec_nd, new_cmp_nd);
        }
    } else {
        new_cmp_nd = new_el->second;
        children_diagnose_old_spec_based_fate(old_spec_nd, new_el->second, top_level);
        record_tax_edits_for_match(old_spec_nd, new_cmp_nd);
    }
}

void TaxonomyDiffer::record_tax_edits_for_match(const RTRichTaxNode * old_nd,
                                                const RTRichTaxNode * new_nd) {
    auto tax_id = old_nd->get_ott_id();
    if (tax_id != new_nd->get_ott_id()) {
        auto & edit = new_alpha_edit();
        edit.operation = AlphaEditOp::CHANGED_ID;
        edit.first_id = tax_id;
        edit.second_id = new_nd->get_ott_id();
    }
    record_syn_diffs(old_nd, new_nd);
    if (old_nd->get_name() != new_nd->get_name()) {
        const auto old_name = old_nd->get_name();
        auto & edit = new_alpha_edit();
        edit.operation = AlphaEditOp::CHANGED_NAME;
        edit.first_str = old_name;
        edit.second_str = new_nd->get_name();
        edit.first_id = tax_id;
    }
    auto & old_data = old_nd->get_data();
    auto & new_data = new_nd->get_data();
    if (old_data.rank != new_data.rank) {
        auto & edit = new_alpha_edit();
        edit.operation = AlphaEditOp::CHANGED_RANK;
        edit.first_rank = old_data.rank;
        edit.second_rank = new_data.rank;
        edit.first_id = tax_id;
    }
    if (old_data.flags != new_data.flags) {
        auto & edit = new_alpha_edit();
        edit.operation = AlphaEditOp::CHANGED_FLAGS;
        edit.first_flags = old_data.flags;
        edit.second_flags = new_data.flags;
        edit.first_id = tax_id;
    }   
}


bool diff_from_taxonomies(std::ostream & out,
                          TaxonomyDiffMaker & old_tax,
                          TaxonomyDiffMaker & new_tax) {
    TaxonomyDiffer tax_dif(old_tax, new_tax);
    tax_dif.write(out);    
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

