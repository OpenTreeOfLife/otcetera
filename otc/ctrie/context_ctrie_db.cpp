#include "otc/ctrie/context_ctrie_db.h"
#include "otc/ctrie/ctrie_db.h"
#include "otc/tnrs/context.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"

using std::set;
using std::string;
using std::vector;

namespace otc {

ContextAwareCTrieBasedDB::ContextAwareCTrieBasedDB(const Context &context_arg,
                                                   const RichTaxonomy &taxonomy)
    :context(context_arg) {
    Context::init_nom_codes_to_traversal(taxonomy);
    if (context_arg.name_matcher != nullptr) {
        return; // already initialized
    }
    const auto & rich_tax_tree = taxonomy.get_tax_tree();
    const auto & rt_data = rich_tax_tree.get_data();
    std::set<std::string> all_names;
    auto insert_hint = all_names.begin();
    for (auto& [name, node] : rt_data.name_to_node)
    {
        // node could be nullptr here if this is a homonym, see note in taxonomy.h
        if (node)
        {
            auto nn = normalize_query(name);
            match_name_to_taxon[nn].push_back(const_rich_taxon_and_syn_ptr{node, nullptr});
            insert_hint = all_names.insert(insert_hint, nn);
        }
    }
    insert_hint = all_names.begin();
    for (auto& [name, nodes] : rt_data.homonym_to_node) {
        auto nn = normalize_query(name);
        for (auto hnp : nodes) {
            assert(hnp);
            match_name_to_taxon[nn].push_back(const_rich_taxon_and_syn_ptr{hnp, nullptr});
        }
        insert_hint = all_names.insert(insert_hint, nn);
    }
    // filtered
    insert_hint = all_names.begin();
    for (auto& [name, record] : rt_data.name_to_record) {
        auto nn = normalize_query(name);
        assert(record);
        match_name_to_taxon[nn].push_back(const_rich_taxon_and_syn_ptr{nullptr, (const void *)record});
        insert_hint = all_names.insert(insert_hint, nn);
    }
    insert_hint = all_names.begin();
    for (auto& [name, records] : rt_data.homonym_to_record) {
        auto nn = normalize_query(name);
        for (auto hrp : records) {
            assert(hrp);
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
    throw OTCError() << "partitioning by context is not implemented yet...";
}


std::set<FuzzyQueryResult, SortQueryResByNearness> ContextAwareCTrieBasedDB::fuzzy_query(const std::string & query_str) const {
    std::set<FuzzyQueryResult, SortQueryResByNearness> sorted;
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

std::set<FuzzyQueryResult, SortQueryResByNearness> ContextAwareCTrieBasedDB::exact_query(const std::string & query_str) const {
    std::set<FuzzyQueryResult, SortQueryResByNearness> sorted;
    if (context.name_matcher != nullptr) {
        sorted = context.name_matcher->exact_query(query_str);
    }
    for (auto c :children) {
        if (c->context.name_matcher) {
            auto csorted = c->context.name_matcher->exact_query(query_str);
            sorted.insert(std::begin(csorted), std::end(csorted));
        }
    }
    return sorted;
}

using vec_fqr_w_t = std::vector<FuzzyQueryResultWithTaxon>;
vec_fqr_w_t ContextAwareCTrieBasedDB::to_taxa(const set<FuzzyQueryResult, SortQueryResByNearness>& sorted,
                                              const RTRichTaxNode * context_root,
                                              const RichTaxonomy & /*taxonomy*/, 
                                              bool include_suppressed) const {
    LOG(DEBUG) << "to_taxa(context_id = " << context_root->get_ott_id() << ", ... , included_suppressed ="  << include_suppressed << ")";
    vec_fqr_w_t results;
    const auto & tax_data = context_root->get_data();
    const auto filter_trav_enter = tax_data.trav_enter;
    const auto filter_trav_exit = tax_data.trav_exit;

    if (sorted.empty()) {
        LOG(DEBUG) << "no matches";
    }
    for (auto fqr : sorted) {
        const auto & vec_taxon_and_syn_ptrs = match_name_to_taxon.at(fqr.match());
        LOG(DEBUG) << "FuzzyQueryResult(match=\"" << fqr.match() << "\", score = " << fqr.score << ") -> vec size = " << vec_taxon_and_syn_ptrs.size();
        for (auto & tax_and_syn_pair : vec_taxon_and_syn_ptrs) {
            auto tax_ptr = tax_and_syn_pair.first;
            if (tax_ptr == nullptr) {
                LOG(DEBUG) << "matched suppressed and include_suppressed = " << include_suppressed;
                if (include_suppressed) {
                    const TaxonomyRecord * tr = (const TaxonomyRecord *)(tax_and_syn_pair.second);
                    results.push_back(FuzzyQueryResultWithTaxon(fqr, tr));
                }
            } else {
                const auto & res_tax_data = tax_ptr->get_data();
                LOG(DEBUG) << "matched taxon trav = (" << res_tax_data.trav_enter <<  ", " << res_tax_data.trav_exit << "). filter.trav = (" << filter_trav_enter << ", " << filter_trav_exit << ")";
                if (res_tax_data.trav_exit <= filter_trav_exit && res_tax_data.trav_enter >= filter_trav_enter) {
                    const TaxonomicJuniorSynonym * syn_ptr = (const TaxonomicJuniorSynonym *)(tax_and_syn_pair.second);
                    if (syn_ptr == nullptr) {
                        LOG(DEBUG) << "pushing non-syn";
                        results.push_back(FuzzyQueryResultWithTaxon(fqr, tax_ptr));
                    } else {
                        LOG(DEBUG) << "pushing synonym";
                        results.push_back(FuzzyQueryResultWithTaxon(fqr, tax_ptr,  syn_ptr));
                    }
                }
            }
        }
    }
    return results;
}

vec_fqr_w_t ContextAwareCTrieBasedDB::fuzzy_query_to_taxa(const std::string & query_str,
                                                          const RTRichTaxNode * context_root,
                                                          const RichTaxonomy & taxonomy,
                                                          bool include_suppressed) const {
    LOG(DEBUG) << "fuzzy_query_to_taxa(" << query_str << ", context_id = " << context_root->get_ott_id() << ", ... , included_suppressed ="  << include_suppressed << ")";
    return to_taxa(fuzzy_query(query_str),  context_root, taxonomy, include_suppressed);
}

} // namespace otc
