#ifndef OTC_WS_NEXSON_NEXSON_H
#define OTC_WS_NEXSON_NEXSON_H
#endif

#include <string>
#include <memory>
#include <map>
#include "otc/tree.h"
#include "ws/tolwsadaptors.h"
#include "json.hpp"

namespace otc
{
    nlohmann::json get_phylesystem_study(const std::string& study_id);
    std::pair<nlohmann::json,nlohmann::json> extract_tree_nexson(const nlohmann::json& nexson, const std::string& treeid);

// https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/org/opentreeoflife/taxa/Nexson.java#120
template<typename T>
std::unique_ptr<T> treeson_get_tree(const nlohmann::json& tree, const nlohmann::json& otus, bool extract_ingroup)
{
    using std::string;

    // 1. Create objects for nodes
    auto nodes = tree["nodeById"];
    std::map<string, typename T::node_type*> node_ptrs;
    typename T::node_type* root = nullptr;
    for(auto x = nodes.begin(); x != nodes.end(); x++)
    {
	auto node = new typename T::node_type(nullptr);
	if (x.value().count("@root"))
	    root = node;

	// set name
	node->set_name(x.key());

	// set ottid
	if (x.value().count("@otu"))
	{
	    string otuid = x.value()["@otu"];
	    auto otu = otus[otuid];
	    if (otu.count("^ot:ottId"))
	    {
		OttId id = otu["^ot:ottId"];
		node->set_ott_id(id);
	    }
	}

	node_ptrs.insert({x.key(), node});
    }

    // 2. Connect source node to target node for each edge
    auto edges = tree["edgeBySourceId"];
    for(auto x = edges.begin(); x != edges.end(); x++)
    {
	const string sourceid = x.key();
	auto source = node_ptrs.at(sourceid);
	// x.value is an object where the keys are edge names and the values are @length
	for(auto& edge: x.value())
	{
	    const string targetid = edge["@target"];
	    auto target = node_ptrs.at(targetid);
	    source->add_child(target);
	}
    }

    // 3. Find root
    auto root_node_id = lookup(tree, "^ot:rootNodeId");
    if (root_node_id)
    {
//	LOG(WARNING)<<"root node = "<<*root_node_id;
	auto root2 = node_ptrs[*root_node_id];
	if (root and root != root2)
	    throw OTCError()<<"Node property '@root' disagrees with tree property '^ot:rootNodeId'!";
	if (not root)
	    root = root2;
    }

    // 4. Make tree from the root (ensure we don't leak non-ingroup nodes)
    auto whole_tree = std::make_unique<T>(root);
    auto ingroup_node_id = lookup(tree, "^ot:inGroupClade");

    // 5. Return ingroup
    if (extract_ingroup and ingroup_node_id)
    {
	auto ingroup_node = node_ptrs[*ingroup_node_id];
	if (not ingroup_node)
	{
	    LOG(WARNING)<<"ingroup_node = NULL  ingroup_id = "<<ingroup_node_id<<"  *ingroup_id = "<<*ingroup_node_id;
	    return whole_tree;
	}
	if (ingroup_node->get_parent()) ingroup_node->detach_this_node();
	return std::make_unique<T>(ingroup_node);
    }
    // 6. Return whole tree
    else
	return whole_tree;
}

// https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/org/opentreeoflife/server/Services.java#L266
template<typename T>
std::unique_ptr<T> get_source_tree(const std::string& study_id, const std::string& tree_id, bool extract_ingroup)
{
    LOG(WARNING)<<"getting phylesystem tree '"<<study_id<<"' '"<<tree_id<<"'";

    // if the study is not found, I think this throws an exception...
    auto study = get_phylesystem_study(study_id);

//    LOG(WARNING)<<"study = '"<<study.dump(1)<<"'";

    auto tree_and_otus = extract_tree_nexson(study, tree_id);

    auto tree = treeson_get_tree<T>(tree_and_otus.first, tree_and_otus.second, extract_ingroup);
    tree->set_name(study_id+"@"+tree_id);

    return tree;
}

// https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/org/opentreeoflife/conflict/ConflictAnalysis.java
// HTTP requests: https://github.com/Corvusoft/restbed/blob/master/example/https_client/source/verify_none.cpp


// https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/org/opentreeoflife/server/Services.java#L242, specToTree(,)
template<typename T>
    std::unique_ptr<T> get_phylesystem_tree(const std::string& study_tree)
{
    using std::list;
    using std::string;
    using std::vector;

    // 1. Determine if delimiter is '#' or '@'
    char delim = '#';
    if (study_tree.find('@') != std::string::npos)
	delim = '@';

    // 2. Decompose study_tree into study@tree
    list<string> parts_list = split_string(study_tree, delim);
    vector<string> parts;
    parts.insert(parts.end(), parts_list.begin(), parts_list.end());

    if (parts.size() != 2)
	throw OTCBadRequest()<<"Expected study@tree or study@tree, but got '"<<study_tree<<"'";

    // 3. Get the study tree from phylesystem
    return get_source_tree<T>(parts[0], parts[1], true);
}

} // namespace otc
