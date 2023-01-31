#include "names.h"

#include "build.h"
#include "oracle.h"
#include "otc/tree_operations.h"

#include <chrono>
#include <map>
#include <utility>
#include <iomanip>

using namespace otc;

using std::vector;
using std::set;
using std::map;
using std::string;

using namespace otc;

Tree_t::node_type* find_mrca_of_desids(const set<OttId>& ids, const std::unordered_map<OttId, Tree_t::node_type*>& summaryOttIdToNode) {
    int first = *ids.begin();
    auto node = summaryOttIdToNode.at(first);
    while( not is_subset(ids, node->get_data().des_ids) )
        node = node->get_parent();
    return node;
}

bool is_ancestor_of(const Tree_t::node_type* n1, const Tree_t::node_type* n2) {
    // make sure the depth fields are initialized
    assert(n1 == n2 or depth(n1) != 0 or depth(n2) != 0);
    if (depth(n2) > depth(n1)) {
        while (depth(n2) != depth(n1)) {
            n2 = n2->get_parent();
        }
        return (n2 == n1);
    } else {
        return false;
    }
}


const Tree_t::node_type* find_unique_maximum(const vector<const Tree_t::node_type*>& nodes)
{
    for(int i=0;i<nodes.size();i++)
    {
        bool is_ancestor = true;
        for(int j=0;j<nodes.size() and is_ancestor;j++) {
            if (j==i) {
                continue;
            }
            if (not is_ancestor_of(nodes[i],nodes[j])) {
                is_ancestor = false;
            }
        }
        if (is_ancestor) {
            return nodes[i];
        }
    }
    return nullptr;
}

Tree_t::node_type* add_monotypic_parent(Tree_t& tree, Tree_t::node_type* nd) {
    if (nd->get_parent()) {
        auto p = nd->get_parent();
        auto monotypic = tree.create_child(p);
        nd->detach_this_node();
        monotypic->add_child(nd);
        return monotypic;
    } else {
        auto monotypic = tree.create_root();
        monotypic->add_child(nd);
        return monotypic;
    }
}

const Tree_t::node_type* select_canonical_ottid(const vector<const Tree_t::node_type*>& names) {
    // We should only have to make this choice if there are at least 2 names to choose from.
    assert(names.size() >= 2);
    // Do something more intelligent here - perhaps prefer non-incertae-sedis taxa, and then choose lowest ottid.
    return names.front();
}

void register_ottid_equivalences(const Tree_t::node_type* canonical, const vector<const Tree_t::node_type*>& names) {
    // First pass - actually we should write on a JSON file.
    std::cerr << canonical->get_name() << " (canonical): equivalent to ";
    for(auto name: names)
        std::cerr << name->get_name() << " ";
    std::cerr << "\n";
}

void add_root_and_tip_names(Tree_t& summary, Tree_t& taxonomy) {
    // name root
    summary.get_root()->set_name(taxonomy.get_root()->get_name());
    if (taxonomy.get_root()->has_ott_id()) {
        summary.get_root()->set_ott_id(taxonomy.get_root()->get_ott_id());
    }
    // name tips
    auto summaryOttIdToNode = get_ottid_to_node_map(summary);
    for(auto nd: iter_leaf(taxonomy)) {
        auto id = nd->get_ott_id();
        auto nd2 = summaryOttIdToNode.at(id);
        nd2->set_name( nd->get_name());
    }
}

//  Given a list of taxon nodes that are known NOT to conflict with the summary tree,
//   * map each name to the MRCA of its include group.
//   * copy the (i) node name and (ii) ottid to the MRCA on the summary tree.
//
//  If multiple names get assigned to a single node, try to handle this by
//    making a monotypic parent assigning the rootmost name to that.
//  Otherwise choose a name arbitrarily.
void add_names(Tree_t& summary, const vector<const Tree_t::node_type*>& compatible_taxa) {
    auto summaryOttIdToNode = get_ottid_to_node_map(summary);

    // 1. Set the des_ids for the summary
    clear_and_fill_des_ids(summary);
    // 2. Place each taxon N at the MRCA of its include group.
    map<Tree_t::node_type*, vector<const Tree_t::node_type*>> name_groups;
    for(auto n2: compatible_taxa) {
        auto mrca = find_mrca_of_desids(n2->get_data().des_ids, summaryOttIdToNode);
        if (not name_groups.count(mrca)) {
            name_groups[mrca] = {};
        }
        name_groups[mrca].push_back(n2);

        // Any extra desids are here because an incertae sedis taxon was placed inside this node,
        // or inside a child.
    }
    // 3. Handle each summary 
    for(auto& name_group: name_groups) {
        auto summary_node = name_group.first;
        auto & names = name_group.second;
        // 3.1. As long as there is a unique root-most name, put that name in a monotypic parent.
        // This can occur when a node has two children, and one of them is an incertae sedis taxon that is moved more tip-ward.
        while (auto max = find_unique_maximum(names)) {
            if (names.size() == 1) {
                set_name_and_maybe_ott_id(*max, *summary_node);
            } else {
                auto p = add_monotypic_parent(summary, summary_node);
                set_name_and_maybe_ott_id(*max, *p);
                p->get_data().des_ids = p->get_first_child()->get_data().des_ids;
            }
            // Move the "removed" elements to the end and the erase them.  Weird.
            names.erase(std::remove(names.begin(), names.end(), max), names.end());
        }
        // 3.2. Select a canonical name from the remaining names.
        if (not names.empty()) {
            // Select a specific ottid as the canonical name for this summary node
            auto canonical = select_canonical_ottid(names);
            summary_node->set_name(canonical->get_name());

            // Write out the equivalence of the remaining ottids to the canonical ottid
            names.erase(std::remove(names.begin(), names.end(), canonical), names.end());
            register_ottid_equivalences(canonical, names);
        }
    }
}

