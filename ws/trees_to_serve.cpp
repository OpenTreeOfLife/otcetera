#include "trees_to_serve.h"

namespace otc
{
    
using std::map;
using std::string;

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
const std::string * TreesToServe::get_stored_string(const std::string & k) {
    const std::string * v = stored_strings[k];
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
using ReadableTaxonomy = std::pair<const RichTaxonomy &,
				   std::unique_ptr<ReadMutexWrapper> >;

using WritableTaxonomy = std::pair<RichTaxonomy &,
				   std::unique_ptr<WriteMutexWrapper> >;

ReadableTaxonomy TreesToServe::get_readable_taxonomy() const {
    assert(taxonomy_ptr != nullptr);
    return {*taxonomy_ptr,
	    std::move(std::make_unique<ReadMutexWrapper>(taxonomy_thread_safety))};
}

WritableTaxonomy TreesToServe::get_writable_taxonomy() {
    assert(taxonomy_ptr != nullptr);
    return {*taxonomy_ptr,
	    std::move(std::make_unique<WriteMutexWrapper>(taxonomy_thread_safety))};
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

TreesToServe::SumTreeInitPair TreesToServe::get_new_tree_and_annotations(const std::string & configfilename,
							   const std::string & filename) {
            
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
    std::unique_ptr<SummaryTree_t> nt = first_newick_tree_from_file<SummaryTree_t>(filename, parsingRules);

    index_by_name_or_id(*nt);
    set_traversal_entry_exit_and_num_tips(*nt);
    tree_list.push_back(move(nt));
    annotation_list.emplace(annotation_list.end());
    auto & sta = annotation_list.back();
    sta.suppressed_from_tree = suppressed_id_set;
    return {*(tree_list.back()), sta};
}

void TreesToServe::register_last_tree_and_annotations() {
    const SummaryTreeAnnotation & sta = annotation_list.back();
    const SummaryTree_t & tree = *(tree_list.back());
    default_synth_id = sta.synth_id; // @ TODO need a better system for deciding the default synth ID.
    id_to_tree[sta.synth_id] = &tree;
    id_to_annotations[sta.synth_id] = &sta;
}

void TreesToServe::free_last_tree_and_annotations() {
    tree_list.back()->clear();
    tree_list.pop_back();
    annotation_list.pop_back();
}

const SummaryTreeAnnotation * TreesToServe::get_annotations(std::string synth_id) const {
    const auto & key = synth_id.empty() ? default_synth_id : synth_id;
    auto mit = id_to_annotations.find(key);
    return mit == id_to_annotations.end() ? nullptr : mit->second;
}

const SummaryTree_t * TreesToServe::get_summary_tree(std::string synth_id) const {
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

    std::map<src_node_id, std::uint32_t> tmpm;
    std::swap(lookup_for_node_ids_while_registering_trees, tmpm);
}


}
