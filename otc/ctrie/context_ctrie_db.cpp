#include "otc/ctrie/context_ctrie_db.h"
#include "otc/ctrie/ctrie_db.h"
#include "otc/tnrs/context.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"

using std::set;
using std::string;
using std::vector;
using std::optional;

namespace otc {

ContextAwareCTrieBasedDB::ContextAwareCTrieBasedDB(const Context &context_arg,
                                                   const RichTaxonomy &taxonomy)
    :context(context_arg) {
    Context::init_nom_codes_boundaries(taxonomy);
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
    for (auto& [name, nodes] : rt_data.homonym_to_nodes) {
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
    Context::init_nom_codes_boundaries(taxonomy);
    throw OTCError() << "partitioning by context is not implemented yet...";
}


std::set<FuzzyQueryResult, SortQueryResByNearness> ContextAwareCTrieBasedDB::fuzzy_query(const std::string & query_str) const {
    std::set<FuzzyQueryResult, SortQueryResByNearness> sorted;
    if (context.name_matcher != nullptr) {
        sorted = context.name_matcher->fuzzy_query(query_str);
    }
    return sorted;
}

optional<string> ContextAwareCTrieBasedDB::exact_query(const std::string & query_str) const
{
    auto nquery = normalize_query(query_str);
    if (match_name_to_taxon.count(nquery))
        return nquery;
    else
        return {};
}

vector<string> ContextAwareCTrieBasedDB::prefix_query(const std::string & query_str) const
{
    auto nquery = normalize_query(query_str);

    if (nquery.size() < 3) return {};

    return trie.prefix_query(nquery);
}

using vec_fqr_w_t = std::vector<FuzzyQueryResultWithTaxon>;
vec_fqr_w_t ContextAwareCTrieBasedDB::to_taxa(const set<FuzzyQueryResult, SortQueryResByNearness>& sorted,
                                              const RTRichTaxNode * context_root,
                                              const RichTaxonomy & /*taxonomy*/, 
                                              bool include_suppressed) const {
    LOG(DEBUG) << "to_taxa(context_id = " << context_root->get_ott_id() << ", ... , included_suppressed ="  << include_suppressed << ")";
    vec_fqr_w_t results;
    const auto & tax_data = context_root->get_data();

    if (sorted.empty()) {
        LOG(DEBUG) << "no matches";
    }
    for (auto fqr : sorted) {
        const auto & vec_taxon_and_syn_ptrs = match_name_to_taxon.at(fqr.match());
//        LOG(DEBUG) << "FuzzyQueryResult(match=\"" << fqr.match() << "\", score = " << fqr.score << ") -> vec size = " << vec_taxon_and_syn_ptrs.size();
        for (auto & [tax_ptr, tax_thing] : vec_taxon_and_syn_ptrs)
        {
            if (tax_ptr == nullptr) {
                LOG(DEBUG) << "matched suppressed and include_suppressed = " << include_suppressed;
                if (include_suppressed) {
                    const TaxonomyRecord * tr = (const TaxonomyRecord *)(tax_thing);
                    results.push_back(FuzzyQueryResultWithTaxon(fqr, tr));
                }
            } else {
                const auto & res_tax_data = tax_ptr->get_data();

                if (is_ancestor_of_using_depth(context_root, tax_ptr))
                {
                    const TaxonomicJuniorSynonym * syn_ptr = (const TaxonomicJuniorSynonym *)(tax_thing);
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

vector<TaxonResult>
ContextAwareCTrieBasedDB::to_taxa(const optional<string>& n_query,
                                  const RTRichTaxNode * context_root,
                                  const RichTaxonomy & /*taxonomy*/, 
                                  bool include_suppressed) const
{
    if (not n_query)
    {
        LOG(DEBUG) << "no matches";
        return {};
    }

    vector<TaxonResult> results;
    const auto & tax_data = context_root->get_data();

    const auto & vec_taxon_and_syn_ptrs = match_name_to_taxon.at(*n_query);
    LOG(DEBUG) << "exact_query(match=\"" << *n_query << ") -> vec size = " << vec_taxon_and_syn_ptrs.size();
    for (auto & [tax_ptr, rec_or_syn_ptr] : vec_taxon_and_syn_ptrs)
    {
        if (tax_ptr == nullptr)
        {
            LOG(DEBUG) << "matched suppressed and include_suppressed = " << include_suppressed;
            if (include_suppressed)
            {
                const TaxonomyRecord * tr = (const TaxonomyRecord *) rec_or_syn_ptr;
                results.push_back(TaxonResult(tr));
            }
        }
        else
        {
            const auto & res_tax_data = tax_ptr->get_data();
            if (is_ancestor_of_using_depth(context_root, tax_ptr))
            {
                const TaxonomicJuniorSynonym * syn_ptr = (const TaxonomicJuniorSynonym *) rec_or_syn_ptr;
                if (syn_ptr == nullptr)
                {
                    LOG(DEBUG) << "pushing non-syn";
                    results.push_back(TaxonResult(tax_ptr));
                }
                else
                {
                    LOG(DEBUG) << "pushing synonym";
                    results.push_back(TaxonResult(tax_ptr,  syn_ptr));
                }
            }
        }
    }
    return results;
}

vector<TaxonResult>
ContextAwareCTrieBasedDB::to_taxa(const vector<string>& n_queries,
                                  const RTRichTaxNode * context_root,
                                  const RichTaxonomy & /*taxonomy*/, 
                                  bool include_suppressed) const
{
    if (n_queries.empty())
    {
        LOG(DEBUG) << "no matches";
        return {};
    }

    vector<TaxonResult> results;

    const auto & tax_data = context_root->get_data();

    for(auto& n_query: n_queries)
    {
        const auto & vec_taxon_and_syn_ptrs = match_name_to_taxon.at(n_query);
        LOG(DEBUG) << "prefix_query(match=\"" << n_query << ") -> vec size = " << vec_taxon_and_syn_ptrs.size();
        for (auto & [tax_ptr, rec_or_syn_ptr] : vec_taxon_and_syn_ptrs)
        {
            if (tax_ptr == nullptr)
            {
                LOG(DEBUG) << "matched suppressed and include_suppressed = " << include_suppressed;
                if (include_suppressed)
                {
                    const TaxonomyRecord * tr = (const TaxonomyRecord *) rec_or_syn_ptr;
                    results.push_back(TaxonResult(tr));
                }
            }
            else
            {
                const auto & res_tax_data = tax_ptr->get_data();
                if (is_ancestor_of_using_depth(context_root, tax_ptr))
                {
                    const TaxonomicJuniorSynonym * syn_ptr = (const TaxonomicJuniorSynonym *) rec_or_syn_ptr;
                    if (syn_ptr == nullptr)
                    {
                        LOG(DEBUG) << "pushing non-syn";
                        results.push_back(TaxonResult(tax_ptr));
                    }
                    else
                    {
                        LOG(DEBUG) << "pushing synonym";
                        results.push_back(TaxonResult(tax_ptr,  syn_ptr));
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

void ContextAwareCTrieBasedDB::add_key(const std::string& s, OttId id, const RichTaxonomy& taxonomy)
{
    auto node = taxonomy.included_taxon_from_id(id);
    if (not node)
	throw OTCError()<<"add_key: id "<<id<<" not found in taxonomy";

    auto nn = normalize_query(s);

    match_name_to_taxon[nn].push_back({node,nullptr});

    trie.add_key(s);
}

} // namespace otc
