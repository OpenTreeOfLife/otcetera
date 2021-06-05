#include "otc/ws/trees_to_serve.h"
#include <regex>

namespace otc
{
    
using std::set;
using std::map;
using std::vector;
using std::string;
using std::pair;
using std::unique_ptr;    

TreesToServe::TreesToServe()
    :taxonomy_thread_safety("taxonomy")
{ }

const src_node_id & TreesToServe::decode_study_node_id_index(std::uint32_t sni_ind) const {
    return src_node_id_storer.at(sni_ind);
}

std::uint32_t TreesToServe::get_source_node_id_index(src_node_id sni) {
    assert(!finalized); // should only be called while registering
    auto it = lookup_for_node_ids_while_registering_trees.find(sni);
    std::uint32_t r;
    if (it == lookup_for_node_ids_while_registering_trees.end()) {
        r = src_node_id_storer.size();
        lookup_for_node_ids_while_registering_trees[sni] = r;
        src_node_id_storer.push_back(sni);
    } else {
        r = it->second;
    }
    return r;
}

const string * TreesToServe::get_stored_string(const string & k) {
    const string * v = stored_strings[k];
    if (v == nullptr) {
        stored_strings_list.push_back(k);
        v = &(stored_strings_list.back());
        stored_strings[k] = v;
    }
    return v;
}

void TreesToServe::set_taxonomy(RichTaxonomy &taxonomy) {
    assert(taxonomy_ptr == nullptr);
    taxonomy_ptr = &taxonomy;
    taxonomy_tree = &(taxonomy.get_tax_tree());
}


ReadableTaxonomy TreesToServe::get_readable_taxonomy() const {
    assert(taxonomy_ptr != nullptr);
    return {*taxonomy_ptr,
            std::make_unique<ReadMutexWrapper>(taxonomy_thread_safety)};
}

WritableTaxonomy TreesToServe::get_writable_taxonomy() {
    assert(taxonomy_ptr != nullptr);
    return {*taxonomy_ptr,
            std::make_unique<WriteMutexWrapper>(taxonomy_thread_safety)};
}

void TreesToServe::fill_ott_id_set(const std::bitset<32> & flags,
                                   OttIdSet & ott_id_set,
                                   OttIdSet & suppressed_from_tree) {
    ott_id_set.clear();
    for (const auto nd : iter_node_const(*taxonomy_tree)) {
        const auto & tax_record_flags = nd->get_data().get_flags();
        auto intersection = flags & tax_record_flags;
        const auto ott_id = nd->get_ott_id();
        if (!intersection.any()) {
            ott_id_set.insert(ott_id);
        } else {
            suppressed_from_tree.insert(ott_id);
        }
    }
}

TreesToServe::SumTreeInitPair TreesToServe::get_new_tree_and_annotations(const string & configfilename,
                                                           const string & filename) {
            
    OttIdSet ott_id_set, suppressed_id_set;
    auto cleaning_flags = cleaning_flags_from_config_file(configfilename);
    fill_ott_id_set(cleaning_flags, ott_id_set, suppressed_id_set);

    assert(taxonomy_ptr != nullptr);

    // Load tree from file
    ParsingRules parsingRules;
    parsingRules.ott_id_validator = &ott_id_set;
    parsingRules.include_internal_nodes_in_des_id_sets = true;
    parsingRules.set_ott_idForInternals = true;
    parsingRules.require_ott_ids = true;
    parsingRules.set_ott_ids = true;
    unique_ptr<SummaryTree_t> nt = first_newick_tree_from_file<SummaryTree_t>(filename, parsingRules);

    index_by_name_or_id(*nt);
    set_traversal_entry_exit_and_num_tips(*nt);
    tree_list.push_back(move(nt));
    annotation_list.emplace(annotation_list.end());
    auto & sta = annotation_list.back();
    sta.suppressed_from_tree = suppressed_id_set;
    return {*(tree_list.back()), sta};
}

vector<int> synth_id_to_version(const string& id) {
    static std::regex VERSION_NUMBER_RE ("(\\d+(\\.\\d+))");
    std::smatch m;
    if (std::regex_search(id, m, VERSION_NUMBER_RE)) {
        vector<int> x;
        for(auto & s: split_string(m[1], '.')) {
            x.push_back(std::stoi(s));
        }
        return x;
    }
    throw OTCError()<<"Synth version '"<<id<<"' has no embedded version number";
}

int compare_versions(const vector<int>& v1, const vector<int>& v2) {
    for(unsigned i = 0; i < std::min(v1.size(),v2.size()); i++) {
        if (v1[i] < v2[i]) {
            return -1;
        }
        if (v1[i] > v2[i]) {
            return 1;
        }
    }
    if (v1.size() < v2.size()) {
        return -1;
    }
    return (v1.size() > v2.size() ? 1 : 0);
}


int compare_synth_ids(const string& id1, const string& id2) {
    auto v1 = synth_id_to_version(id1);
    auto v2 = synth_id_to_version(id2);
    return compare_versions(v1,v2);
}

void TreesToServe::register_last_tree_and_annotations() {
    const SummaryTreeAnnotation & sta = annotation_list.back();
    const SummaryTree_t & tree = *(tree_list.back());
    id_to_tree[sta.synth_id] = &tree;
    id_to_annotations[sta.synth_id] = &sta;
    if (default_synth_id.empty() or compare_synth_ids(default_synth_id,sta.synth_id) < 0) {
        default_synth_id = sta.synth_id;
    }
}

void TreesToServe::free_last_tree_and_annotations() {
    tree_list.back()->clear();
    tree_list.pop_back();
    annotation_list.pop_back();
}

string TreesToServe::get_default_tree() const
{
    return default_synth_id;
}

set<string> TreesToServe::get_available_trees() const
{
    set<string> synth_ids;
    for(auto& x: id_to_tree) {
        synth_ids.insert(x.first);
    }
    return synth_ids;
}

const SummaryTreeAnnotation * TreesToServe::get_annotations(string synth_id) const {
    const auto & key = synth_id.empty() ? default_synth_id : synth_id;
    LOG(DEBUG) << " key = \"" << key << "\" synth_id = \"" << synth_id << "\"";
    auto mit = id_to_annotations.find(key);
    if (mit == id_to_annotations.end()) {
        return nullptr;
    }
    return mit->second;
}

const SummaryTree_t * TreesToServe::get_summary_tree(string synth_id) const {
    const auto & key = synth_id.empty() ? default_synth_id : synth_id;
    auto mit = id_to_tree.find(key);
    return mit == id_to_tree.end() ? nullptr : mit->second;
}

std::size_t TreesToServe::get_num_trees() const {
    return id_to_tree.size();
}

void TreesToServe::final_tree_added() {
    finalized = true;
    if (get_num_trees() == 1) {
        const auto & annot = annotation_list.back();
        const auto & sft = annot.suppressed_from_tree;
        assert(taxonomy_ptr != nullptr);
        // TODO: make taxonomy_ptr non-const or make the relevant field mutable. Probably the former.
        auto nct = const_cast<RichTaxonomy *>(taxonomy_ptr);
        nct->set_ids_suppressed_from_summary_tree_alias(&sft);
    }

    map<src_node_id, std::uint32_t> tmpm;
    std::swap(lookup_for_node_ids_while_registering_trees, tmpm);
}


}
