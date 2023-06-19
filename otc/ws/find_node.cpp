#include "find_node.h"
#include <regex>
#include "otc/ws/node_namer_supported_by_stasher.h"

using std::vector;
using std::map;
using std::pair;
using std::tuple;
using std::set;
using std::string;
using std::string_view;
using std::optional;
using std::ostringstream;
using std::unique_ptr;

using namespace boost::property_tree;
using json=nlohmann::json;


using otc::OttId;

namespace otc {

const std::regex ott_id_pattern("^ott(\\d+)$");

optional<OttId> is_ott_id(const string& node_id)
{
    std::smatch matches;
    if (std::regex_match(node_id, matches, ott_id_pattern))
    {
        long raw_ott_id = long_ott_id_from_name(node_id);
        if (raw_ott_id >= 0)
            return to_OttId(raw_ott_id);
    }
    return {};
}

const std::regex mrca_id_pattern("^mrca(ott\\d+)(ott\\d+)$");

TaxonToSynth find_node_by_valid_ottid(const SummaryTree_t & tree, OttId id, const string& node_id)
{
    const auto & tree_data = tree.get_data();

    // 1. Try to find the OTT ID in the summary tree.
    auto i2nit = tree_data.id_to_node.find(id);
    if (i2nit != tree_data.id_to_node.end())
        return TaxonFound{i2nit->second};

    // 2. We didn't find a summary tree node for this OTT ID.  Is this taxon listed as broken?
    if (auto bt_it = tree_data.broken_taxa.find(node_id); bt_it != tree_data.broken_taxa.end())
    {
        // if found it listed as broken, we return the MRCA pointer in the first slot of the pair.
        return TaxonBroken{bt_it->second.first};
    }
    else
    {
        LOG(WARNING) << "OTT ID" << id << " (from '"<<node_id<<"') is not in the synth tree, and is not listed as broken.";
        return TaxonPruned{};
    }
}

OTTNameToSynth find_node_by_ottid_str(const SummaryTree_t & tree, const RichTaxonomy& taxonomy, const string & node_id)
{
    const auto & tree_data = tree.get_data();

    std::smatch matches;
    assert(std::regex_match(node_id, matches, ott_id_pattern));

    // Get the OTT ID
    auto ott_id = is_ott_id(node_id);

    // 1. Check that the ID is not too large
    if (not ott_id)
    {
        LOG(WARNING) << "OTT ID from "<<node_id<<"' is too large!";
        return BadID{};
    }

    // 2. Try and forward the ID.
    auto valid_ott_id = taxonomy.get_unforwarded_id(*ott_id);
    if (not valid_ott_id)
    {
        LOG(WARNING) << "OTT ID " << *ott_id << " (from '"<<node_id<<"') is neither a current ID nor a forwarded ID.";
        return InvalidID{*ott_id};
    }
    auto forwarded_from = ott_id;
    if (*ott_id == *valid_ott_id)
        forwarded_from = {};

    // 3. Map the valid ottid to the summary tree.
    return OTTNameToSynth{ValidID{*valid_ott_id, forwarded_from, find_node_by_valid_ottid(tree, *valid_ott_id, node_id)}};
}

MRCANameToSynth find_node_by_mrca_str(const SummaryTree_t & tree, const RichTaxonomy& taxonomy, const string & node_id)
{
    const auto & tree_data = tree.get_data();

    std::smatch matches;
    auto match = std::regex_match(node_id, matches, mrca_id_pattern);
    assert(match);
    assert(matches.size() >= 2);

    // We used to look up the canonical name, which would be faster...
    //        auto n2nit = tree_data.broken_name_to_node.find(node_id);
    //        if (n2nit != tree_data.broken_name_to_node.end()) {
    //            return {n2nit->second};
    //        }

    std::string first_id = matches[1];
    std::string second_id = matches[2];
    auto result1 = find_node_by_ottid_str(tree, taxonomy, first_id);
    auto result2 = find_node_by_ottid_str(tree, taxonomy, second_id);

    const SumTreeNode_t* mrca = nullptr;
    if (result1.node() and result2.node())
        mrca = find_mrca_via_traversal_indices(result1.node(), result2.node());

    return MRCANameToSynth{result1, result2, mrca};
}

NameToSynth find_node_by_id_str(const SummaryTree_t & tree, const RichTaxonomy& taxonomy, const string & node_id)
{
    std::smatch matches;

    if (std::regex_match(node_id, matches, ott_id_pattern))
        return find_node_by_ottid_str(tree, taxonomy, node_id);

    else if (std::regex_match(node_id, matches, mrca_id_pattern))
        return find_node_by_mrca_str(tree, taxonomy, node_id);

    else
        return NoMatchName{};
}

// Possible statuses:
//  "unknown_id"        (The number in the ottid is too big)
//  "unknown_id"        (Deprecated: previously valid ottid, no longer forwarded)
//  "unknown_id"        (For an mrcaottXottY id where anything goes wrong at all)
//  "invalid_ott_id"    (Not in current ott, and not forwarded)
//  "pruned_ott_id"     (In OTT, but pruned from synth)

//  "broken"            (In OTT, not pruned, but broken taxon and fail_broken = true)
string find_node_failure_reason(const NameToSynth& result, bool fail_broken)
{
    assert(not result.node() or (fail_broken and result.broken()));

    if (result.invalid())
        return "invalid_ott_id";
    else if (result.pruned())
        return "pruned_ott_id";
    else if (result.broken())
        return "broken";
    else
        return "unknown_id";
}


NameToSynth find_required_node_by_id_str(const SummaryTree_t & tree, const RichTaxonomy& taxonomy, const string & node_id)
{
    auto result = find_node_by_id_str(tree, taxonomy, node_id);
    if (not result.node())
    {
        string reason = find_node_failure_reason(result);
        throw OTCBadRequest() << "node_id '" << node_id << "' was not found!"<<json{ {"reason", reason} };
    }
    return result;
}

tuple<vector<const SumTreeNode_t*>,json,json>
find_nodes_for_id_strings(const RichTaxonomy& taxonomy, const SummaryTree_t* tree_ptr, const vector<string>& node_ids,
                          bool fail_broken, bool filter_invalid, bool filter_pruned, bool filter_broken)
{
    vector<const SumTreeNode_t *> nodes;
    json unknown;
    json broken = json::object();
    json filtered = json::object();
    optional<string> bad_node_id;
    for (auto node_id : node_ids)
    {
        auto result = find_node_by_id_str(*tree_ptr, taxonomy, node_id);

        if (not result.node() or (result.broken() and (fail_broken or filter_broken)))
        {
            // Possible statuses:
            //  "unknown_id"        (The number in the ottid is too big)
            //  "unknown_id"        (Deprecated: previously valid ottid, no longer forwarded)
            //  "unknown_id"        (For an mrca where anything goes wrong at all"
            //  "invalid_ott_id"    (Not in current ott, and not forwarded)
            //  "pruned_ott_id"     (In OTT, but pruned from synth)
            //  "broken"            (In OTT, not pruned, but broken taxon and fail_broken = true)

            string reason = "unknown_id";
            bool filter = false;

            if (result.invalid())
            {
                reason = "invalid_ott_id";
                filter = filter_invalid;
            }
            else if (result.pruned())
            {
                reason = "pruned_ott_id";
                filter = filter_pruned;
            }
            else if (result.broken())
            {
                reason = "broken";
                filter = filter_broken;
            }

            if (filter)
                filtered[node_id] = reason;
            else
            {
                unknown[node_id] = reason;
                bad_node_id = node_id;
            }
        }

        if (result.broken())
            broken[node_id] = node_id_for_summary_tree_node(*result.node());

        // Current default strategy means that we include MRCAs for broken taxa.
        nodes.push_back(result.node());
    }
    if (unknown.size()) {
        // Also add {"filtered", filtered} here, if we decide to allow filtering.
        // Probability we should report "filtered": {} when nothing is filtered, but could report filtered only if non-empty.
        throw OTCBadRequest()<<"node_id '"<< *bad_node_id << "' was not found!"<<json{ {"unknown", unknown} };
    }
    return {nodes, broken, filtered};
}

}
