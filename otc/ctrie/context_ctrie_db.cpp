#include "otc/ctrie/context_ctrie_db.h"
#include "otc/ctrie/ctrie_db.h"
#include "otc/tnrs/context.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"

namespace otc {

ContextAwareCTrieBasedDB::ContextAwareCTrieBasedDB(const Context &context_arg,
                                                   const RichTaxonomy &taxonomy)
    :context(context_arg) {
    Context::init_nom_codes_to_traversal(taxonomy);
    init_char_maps();
    if (context_arg.name_matcher != nullptr) {
        return; // already initialized
    }
    const auto & rich_tax_tree = taxonomy.get_tax_tree();
    const auto & rt_data = rich_tax_tree.get_data();
    std::set<std::string> all_names;
    auto insert_hint = all_names.begin();
    for (auto const & name2nd : rt_data.name_to_node) {
        auto nn = normalize_query(name2nd.first);
        match_name_to_taxon[nn].push_back(const_rich_taxon_and_syn_ptr{name2nd.second, nullptr});
        insert_hint = all_names.insert(insert_hint, nn);
    }
    insert_hint = all_names.begin();
    for (auto name2ndvec : rt_data.homonym_to_node) {
        auto nn = normalize_query(name2ndvec.first);
        for (auto hnp : name2ndvec.second) {
            match_name_to_taxon[nn].push_back(const_rich_taxon_and_syn_ptr{hnp, nullptr});
        }
        insert_hint = all_names.insert(insert_hint, nn);
    }
    // filtered
    insert_hint = all_names.begin();
    for (auto name2rec : rt_data.name_to_record) {
        auto nn = normalize_query(name2rec.first);
        match_name_to_taxon[nn].push_back(const_rich_taxon_and_syn_ptr{nullptr, (const void *)name2rec.second});
        insert_hint = all_names.insert(insert_hint, nn);
    }
    insert_hint = all_names.begin();
    for (auto name2recvec : rt_data.homonym_to_record) {
        auto nn = normalize_query(name2recvec.first);
        for (auto hrp : name2recvec.second) {
            match_name_to_taxon[nn].push_back(const_rich_taxon_and_syn_ptr{nullptr, (const void *)hrp});
        }
        insert_hint = all_names.insert(insert_hint, nn);
    }
    for (const auto & tjs : taxonomy.get_synonyms_list()) {
        auto nn = normalize_query(tjs.name);
        match_name_to_taxon[nn].push_back(const_rich_taxon_and_syn_ptr{tjs.primary, (const void *)(&tjs)});
        all_names.insert(nn);
    }
    
    trie.initialize(all_names);
    context_arg.name_matcher = &trie;
}


ContextAwareCTrieBasedDB::ContextAwareCTrieBasedDB(const Context &context_arg,
                         const RichTaxonomy & taxonomy,
                         const std::set<std::string_view> &)
    :context(context_arg) {
    Context::init_nom_codes_to_traversal(taxonomy);
    init_char_maps();
    throw OTCError() << "partitioning by context is not implemented yet...";
}


sorted_q_res_set ContextAwareCTrieBasedDB::fuzzy_query(const std::string & query_str) const {
    sorted_q_res_set sorted;
    if (context.name_matcher != nullptr) {
        sorted = context.name_matcher->fuzzy_query(query_str);
    }
    for (auto c :children) {
        if (c->context.name_matcher) {
            auto csorted = c->context.name_matcher->fuzzy_query(query_str);
            sorted.insert(std::begin(csorted), std::end(csorted));
        }
    }
    return sorted;
}

std::optional<FuzzyQueryResult>  ContextAwareCTrieBasedDB::exact_query(const std::string & raw_query, const std::string & norm_query) const {
    if (context.name_matcher != nullptr) {
        auto res = context.name_matcher->exact_query(raw_query, norm_query);
        if (res) {
            return res;
        }
    }
    for (auto c :children) {
        if (c->context.name_matcher) {
            auto res = c->context.name_matcher->exact_query(raw_query, norm_query);
            if (res) {
                return res;
            }
        }
    }
    return {};
}

struct SortQueryResWTaxonByNearness {
    bool operator() (const FuzzyQueryResultWithTaxon & lhs,
                     const FuzzyQueryResultWithTaxon & rhs) const {
        if (lhs.get_score() < rhs.get_score()) {
            return false;
        } else if (rhs.get_score() < lhs.get_score()) {
            return true;
        }
        return lhs.get_matched_name() < rhs.get_matched_name();
    }
};



vec_q_res_w_taxon ContextAwareCTrieBasedDB::tie_to_taxa(const sorted_q_res_set & sorted,
                                                  const std::string & query_str,
                                                  const RTRichTaxNode * context_root,
                                                  const RichTaxonomy & , 
                                                  keep_taxon_pred_t keep,
                                                  const std::string * exact_string) const {
    const auto & tax_data = context_root->get_data();
    const auto filter_trav_enter = tax_data.trav_enter;
    const auto filter_trav_exit = tax_data.trav_exit;
    if (sorted.empty()) {
        LOG(DEBUG) << "no matches";
    }
    const auto qlc = lower_case_version(query_str);
    const auto wc = to_u32string(qlc);
    auto wcp = &wc;
    std::set<FuzzyQueryResultWithTaxon, SortQueryResWTaxonByNearness> sorted_correct_score;
    for (auto fqr : sorted) {
        const auto & vec_taxon_and_syn_ptrs = match_name_to_taxon.at(fqr.match());
        LOG(DEBUG) << "FuzzyQueryResult(match=\"" << fqr.match() << "\", score = " << fqr.score << ") -> vec size = " << vec_taxon_and_syn_ptrs.size();
        for (auto & tax_and_syn_pair : vec_taxon_and_syn_ptrs) {
            auto tax_ptr = tax_and_syn_pair.first;
            if (not keep(tax_ptr)) {
                LOG(DEBUG) << "skipped based on predicate returning false";
                continue;
            }
            if (tax_ptr == nullptr) {
                LOG(DEBUG) << "taxon record match";
                const TaxonomyRecord * tr = (const TaxonomyRecord *)(tax_and_syn_pair.second);
                if (exact_string == nullptr || tr->name == *exact_string) {
                    sorted_correct_score.emplace(FuzzyQueryResultWithTaxon{fqr, tr, wcp});
                }
           } else {
                const auto & res_tax_data = tax_ptr->get_data();
                LOG(DEBUG) << "matched taxon trav = (" << res_tax_data.trav_enter <<  ", " << res_tax_data.trav_exit << "). filter.trav = (" << filter_trav_enter << ", " << filter_trav_exit << ")";
                if (res_tax_data.trav_exit <= filter_trav_exit && res_tax_data.trav_enter >= filter_trav_enter) {
                    const TaxonomicJuniorSynonym * syn_ptr = (const TaxonomicJuniorSynonym *)(tax_and_syn_pair.second);
                    if (syn_ptr == nullptr) {
                        LOG(DEBUG) << "pushing non-syn";
                        if (exact_string == nullptr || tax_ptr->get_name() == *exact_string) {
                            sorted_correct_score.emplace(FuzzyQueryResultWithTaxon{fqr, tax_ptr, wcp});
                        }
                    } else {
                        LOG(DEBUG) << "pushing synonym";
                        if (exact_string == nullptr || syn_ptr->get_name() == *exact_string) {
                            sorted_correct_score.emplace(FuzzyQueryResultWithTaxon{fqr, tax_ptr,  syn_ptr, wcp});
                        }
                    }
                }
            }
        }
    }
    vec_q_res_w_taxon results;
    for (const auto & sr : sorted_correct_score) {
        results.push_back(sr);
    }
    return results;
}
    

} // namespace otc
