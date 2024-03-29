#include <regex>
#include "otc/ws/tolws.h"
#include "otc/ws/tolwsadaptors.h"
#include "otc/ws/trees_to_serve.h"
#include "otc/conflict.h"
#include "otc/ws/nexson/nexson.h"
#include "otc/ws/prune.h"
#include "otc/ws/find_node.h"
#include <optional>
#include <string_view>


using std::vector;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::string_view;
using std::optional;
using std::ostringstream;
using std::unique_ptr;

using namespace boost::property_tree;
using json=nlohmann::json;


namespace otc {

using cnode_type = ConflictTree::node_type;
using snode_type = SummaryTree_t::node_type;



// Get the subtree of T1 connecting the leaves of T1 that are also in T2.
template <typename Tree_Out_t, typename Tree1, typename Tree2>
pair<unique_ptr<Tree_Out_t>,unique_ptr<Tree_Out_t>>
get_induced_trees2(Tree1& T1,
                   std::function<const typename Tree1::node_type*(const typename Tree1::node_type*,const typename Tree1::node_type*)> MRCA_of_pair1,
                   Tree2& T2,
                   std::function<const typename Tree2::node_type*(const typename Tree2::node_type*,const typename Tree2::node_type*)> MRCA_of_pair2)
{
    LOG(DEBUG)<<"T1 = "<<newick_string(T1);
    LOG(WARNING)<<"n_leaves(T1) = "<<n_leaves(T1);
    LOG(WARNING)<<"n_leaves(T2) = "<<n_leaves(T2);

    // 1. First construct the induced tree for T2.

    // 1a. Find the nodes of T2 that correspond to leaves of T1.
    //     Note that some of these nodes could be ancestral to other ones in T2.
    auto T2_nodes_from_T1_leaves = get_induced_nodes(T1,T2);
    LOG(WARNING)<<T2_nodes_from_T1_leaves.size()<<" leaves of T1 are in T2";

    if (T2_nodes_from_T1_leaves.size() < 2)
        throw OTCBadRequest()<<"The two trees have only "<<T2_nodes_from_T1_leaves.size()<<" nodes in common!  At least 2 are required";

    // 1b. Actually construct the induced tree for T2.
    //     It might have fewer leaves than in T2_nodes_from_T1_leaves, if some of the nodes are ancestral to others.
    auto induced_tree2 = get_induced_tree<Tree_Out_t>(T2_nodes_from_T1_leaves, MRCA_of_pair2);
    LOG(WARNING)<<"n_leaves(induced_tree2) = "<<n_leaves(*induced_tree2);
    LOG(DEBUG)<<"induced-tree2a = "<<newick_string(*induced_tree2);

    // 1c. Rename internal nodes of synth/taxonomy to ottXXX instead of taxon name
    for(auto node: iter_post(*induced_tree2)) {
        // ott: Parietobalaena palmeri -> ott3615461
        // synth: "" -> ott494235
        // newick: "_node41 ott3612189" -> ott3612189

        // Its only newick trees from phylesystem that have source node names, so we
        // can safely use the source node name if its present.
        if (auto source_name = get_source_node_name(node->get_name()))
            node->set_name(*source_name);
        else if (node->has_ott_id()) {
            string new_name = "ott"+std::to_string(node->get_ott_id());
            LOG(DEBUG)<<"Renaming "<<node->get_name()<<" to "<<new_name;
            node->set_name(new_name);
        }
    }
    LOG(DEBUG)<<"induced_tree2 = "<<newick_string(*induced_tree2);

    //FIXME - handle cases like (Homo sapiens, Homo) by deleting monotypic nodes at the root.

    // 2. Keep only leaves from T1 that (a) have an ottid in T2 and (b) map to a leaf of induced_tree2
    auto ottid_to_induced_tree2_node = otc::get_ottid_to_const_node_map(*induced_tree2);
    std::vector<const typename Tree1::node_type*> induced_leaves1;
    for(auto leaf: iter_leaf_const(T1))
    {
        if (not leaf->has_ott_id())
        {
            LOG(DEBUG)<<"Dropping tip: no ott id";
            continue;
        }

        auto it = ottid_to_induced_tree2_node.find(leaf->get_ott_id());
        if (it == ottid_to_induced_tree2_node.end()) {
            LOG(DEBUG)<<"Dropping query tip "<<leaf->get_ott_id()<<": not found in induced tree2.";
        } else if (not it->second->is_tip()) {
            LOG(DEBUG)<<"Dropping higher taxon tip "<<leaf->get_ott_id();
        } else {
            induced_leaves1.push_back(leaf);
        }
    }
    LOG(WARNING)<<"keeping "<<induced_leaves1.size()<<" leaves from T1";

    // 3. Construct the induced tree for T1
    auto induced_tree1 = get_induced_tree<Tree_Out_t>(induced_leaves1, MRCA_of_pair1);
    LOG(WARNING) << "n_leaves(induced_tree1) = " << n_leaves(*induced_tree1);
    LOG(DEBUG) << "induced_tree1 = " << newick_string(*induced_tree1);
    assert(n_leaves(*induced_tree1) == n_leaves(*induced_tree2));
    return {std::move(induced_tree1), std::move(induced_tree2)};
}


template <typename N>
inline bool is_fake_tip(const N* n)
{
    return (n->get_out_degree() > 0) and (n->get_first_child()->get_name().empty());
}

// input tree node names should all be of the form node###
// ott        node names should all be of the form ott###
// synth tree node names should all be of the form ott###  or  mrcaott###ott###
template <typename N>
string node_name(const N* node) {
    assert(not node->get_name().empty());
    return node->get_name();
}

int min_depth(const cnode_type* node) {
    while (node and node->get_parent() and node->get_parent()->is_outdegree_one_node()) {
        node = node->get_parent();
    }
    return node->get_data().depth;
}

typedef std::function<optional<string>(const string&)> witness_namer_t;

/*
 * This structure exists so that conflict analyis can record results on it.
 */
struct conflict_stats {
    protected:
    template <typename T>
    void throw_otcerror_if_second_not_true(const T & result, const cnode_type* node1) {
        if (not result.second) {
            throw OTCError() << "key "<<node_name(node1) << " occurs twice in a conflict_stats map!";
        }    
    }
    public:

    /*
     * We store strings here instead of (const cnode_type*) because the nodes
     * may vanish out from underneath us.  But pointers would be nicer.
     */

    witness_namer_t witness_namer;

    map<string, string> supported_by;
    map<string, string> partial_path_of;
    map<string, std::set<pair<string, int>>> conflicts_with;
    map<string, string> resolved_by;
    map<string, string> resolves;
    map<string, string> terminal;

    void add_supported_by(const cnode_type* node2, const cnode_type* node1) {
        auto result = supported_by.insert({node_name(node1), node_name(node2)});
        // We can have the relationship node2 `supported_by` node1 for only one node2.
        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_partial_path_of(const cnode_type* node2, const cnode_type* node1) {
        partial_path_of.insert({node_name(node1), node_name(node2)});
        // We can have the relationship node2 `partial_path_of` node1 for multiple node2s.
        // throw_otcerror_if_second_not_true(result, node1);
    }

    void add_resolved_by(const cnode_type* node2, const cnode_type* node1) {
        auto result = resolved_by.insert({node_name(node1), node_name(node2)});
        // We can have the relationship node2 `resolved_by` node1 for only one node2.
        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_resolves(const cnode_type* node2, const cnode_type* node1) {
        resolves.insert({node_name(node1), node_name(node2)});
        // We can only have the relationship (node2 `resolves` node1) one multiple node2s.
        // throw_otcerror_if_second_not_true(result, node1);
    }

    void add_terminal(const cnode_type* node2, const cnode_type* node1) {
        terminal.insert({node_name(node1), node_name(node2)});
//      We can have the relationship (node2 `terminal` node1) for multiple node2s!
//      Let's keep the FIRST one, since it will be the most tipward-one.
//        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_conflicts_with(const cnode_type* node2, const cnode_type* node1)
    {
        auto name1 = node_name(node1);
        auto name2 = node_name(node2);

        // 1. Get the (possibly empty) set of deepest nodes in tree2 conflicting with node1 in tree1.
        set<pair<string, int>>& nodes = conflicts_with[name1];

        int min_depth_node2 = min_depth(node2);

        // 2. If the newest node (node2) is deeper, then clear the non-minimal-depth nodes.
        if (not nodes.empty() and (min_depth_node2 < nodes.begin()->second)) {
            nodes.clear();
        }
        // 3. If the newest node (node2) is shallower then do nothing
        if (not nodes.empty() and (min_depth_node2 > nodes.begin()->second)) {
            return;
        }
        // 4. If the newest node is not shallower then add it.
        nodes.insert(pair<string, int>({name2, min_depth_node2}));
    }
    json get_json(const ConflictTree& induced_tree1) const;

    conflict_stats(const witness_namer_t& w): witness_namer(w) {}
};

using tnode_type = RTRichTaxNode;

inline int depth(const tnode_type* node) {
    return node->get_data().depth;
}

string extract_node_name_if_present(const string& newick_name) {
    auto node_name = get_source_node_name(newick_name);
    if (node_name) {
        return *node_name;
    }
    return newick_name;
}

optional<string> ott_witness_namer(const string& witness, const RichTaxonomy& Tax)
{
    long raw_ott_id = long_ott_id_from_name(witness);
    if (raw_ott_id >= 0) {
        OttId id = check_ott_id_size(raw_ott_id);
        auto nd = Tax.included_taxon_from_id(id);
        if (nd != nullptr) {
            return nd->get_name();
        }
    }
    return {};
}

optional<pair<string,int>> get_descendant_name(const RichTaxonomy& taxonomy, const snode_type* nd)
{
    // 1. If the node has a name, the use that name.
    if (nd->has_ott_id())
    {
        // FIXME -- should we use ott_witness_namer(node_name(*nd)) instead?
        auto name = taxon_nonuniquename(taxonomy, *nd);
        return {{name, 0}};
    }

    // 2. Otherwise consider names of its child subtrees
    optional<pair<string,int>> name;
    for(auto child: iter_child_const(*nd))
    {
        if (auto child_name = get_descendant_name(taxonomy, child))
        {
            child_name->second++;
            if (not name or child_name->second < name->second)
                name = child_name;
        }
    }
    return name;
}

// Compare to get_descendant_names( ) in tolws.cpp
vector<string> get_descendant_names2(const RichTaxonomy& taxonomy, const snode_type* nd)
{
    vector<pair<string,int>> weighted_names;
    for(auto child: iter_child_const(*nd))
    {
        if (auto n = get_descendant_name(taxonomy, child))
            weighted_names.push_back(*n);
    }

    std::sort(weighted_names.begin(), weighted_names.end(), [](auto& x, auto& y) {return x.second < y.second;});

    vector<string> names;
    for(auto& [name,_]: weighted_names)
        names.push_back(name);
    return names;
}

optional<string> synth_witness_namer(const string& node_name, const SummaryTree_t& summary, const RichTaxonomy& taxonomy)
{
    // 1. If we can name the node using the taxonomy, then do that.
    if (auto wn = ott_witness_namer(node_name, taxonomy))
        return *wn;

    // 2. Look up the node in the summary tree.
    auto node = find_required_node_by_id_str(summary, taxonomy, node_name).node();

    // 3. Find the names of descendants
    auto dnames = get_descendant_names(taxonomy, *node);

    // 4. If there are 0 or 1 name, then we didn't find anything to call this node.
    if (dnames.size() < 2)
        return {};

    // 5. Generate the witness_name
    string name = dnames.front() + "+" + dnames.back();
    if (dnames.size() > 2)
        name += "+...";
    name = "["+name+"]";

    return name;
}

inline json get_node_status(const string& witness, string status, const witness_namer_t& witness_namer) {
    json j;
    j["witness"] = extract_node_name_if_present(witness);
    j["status"] = std::move(status);
    if (auto w = witness_namer(witness))
        j["witness_name"] = *w;
    return j;
}

json get_conflict_node_status(const set<pair<string,int>>& witnesses,
                              string status,
                              const witness_namer_t& witness_namer)
{
    json j_witnesses = json::array();
    json j_witness_names = json::array();
    for(auto& [node_name, _]: witnesses)
    {
        json witness_name;
        if (auto w = witness_namer(node_name))
            witness_name = *w;

        // The witness could be
        // (i) an ottX name (tree2 = ott or synth)
        // (ii) an mrcaottXottY name (tree2 = synth)
        // (iii) a nodeX name (tree2 = study tree)

        // Can we provide a meaningful name for case (ii)?

        j_witnesses.push_back( extract_node_name_if_present(node_name) );
        j_witness_names.push_back( witness_name );
    }

    json j;
    j["status"] = std::move(status);
    j["witness_name"] = j_witness_names;
    j["witness"] = j_witnesses;
    return j;
}

// FIXME: Conflict relations currently refer to nodes by a name string.
//        We assume that this name string:
//         1. Is accessed as node->get_name()
//         2. Contains the ottid, if there is one.
//        The curator application also assumes that
//         3. JSON results reference the nexml node name (e.g. nodeYYY) for phylesystem trees.
//
//        Actual node->get_name() strings look like this:
//         ott:   ottXXX
//         synth: ottXXX or mrcaottWWWottZZZ
//         (internally) phylesystem: nodeYYY
//         newick: stuff_ottXXX for leaves, and stuff otherwise.
//         newick from phylesystem (externally): nodeYYY (for internal nodes) or nodeYYY_ottXXX (for leaves).
//
//        I assume that phylesystem node names are transformed from name -> name_ottXXX for leaf nodes, when translating phylesystem trees to newick.
//
//        Assumptions 1-3 are met for ott, synth, and (internal) phylesystem trees.
//
//        For newick trees from phylesystem, we handle this situation by translating these names from
//        stuff_nodeYYY_ottXXX -> nodeYYY before constructing the json, as below. (10/03/2017 - BDR).
//
//        *ALTERNATIVELY*, we could chop any _ottXXX suffix from the leaf names to get the name used in the JSON result, but this could affect non-phylesystem inputs.
//


// PROBLEM: It is possible to have both y1 conflicts with x (x:conflicts_with y1) and y2 resolves x (x:resolves y2)
// PROBLEM: It is possible to have both y1 supported_by x (x:supported_by y1) and y2 resolves x (x:resolves y2)
// PROBLEM: It is possible to have both x resolves y (x:resolved_by y1) and y2 resolves x (x:resolves y2)
// SOLUTION: DO NOT report when tree2 resolves tree1, only when tree1 resolves tree2.


json conflict_stats::get_json(const ConflictTree& tree1) const
{
    json nodes;

/*
 *  node1: node from tree1.
 *  node2: node from tree2.
 *
 *  NOTE: All the relationships describe the nodes in tree2!
 *        For example, "resolved_by" means that node2 is resolved_by node1.
 *
 *        Despite this fact, the map is keyed by nodes in the FIRST tree.
 *        Leaving out the "tree2 resolves tree1" relationship allows us to divide
 *          nodes in tree1 by how they relate to tree2.
 */ 

//    for(auto& [node1, node2]: resolves) {
//        nodes[extract_node_name_if_present(node1)] = get_node_status(node2, "resolves", witness_namer);
//    }
    for(auto& [node1, node2]: resolved_by) {
        nodes[extract_node_name_if_present(node1)] = get_node_status(node2, "resolved_by", witness_namer);
    }
    for(auto& [node1, node2]: supported_by) {
        nodes[extract_node_name_if_present(node1)] = get_node_status(node2, "supported_by", witness_namer);
    }
    for(auto& [node1, node2]: partial_path_of) {
        nodes[extract_node_name_if_present(node1)] = get_node_status(node2, "partial_path_of", witness_namer);
    }
    for(auto& [node1, node2]: terminal) {
        nodes[extract_node_name_if_present(node1)] = get_node_status(node2, "terminal", witness_namer);
    }
    for(auto& [node1, node2]: conflicts_with) {
        nodes[extract_node_name_if_present(node1)] = get_conflict_node_status(node2, "conflicts_with", witness_namer);
    }
    // For monotypic nodes in the query, copy annotation from child.
    for(auto it: iter_post_const(tree1)) {
        if (it->is_outdegree_one_node()) {
            auto name = extract_node_name_if_present(it->get_name());
            auto child_name = extract_node_name_if_present(it->get_first_child()->get_name());
            assert(name.size());
            assert(child_name.size());
            nodes[name] = nodes.at(child_name);
        }
    }
    return nodes;
}

/*
 * other_tree is (currently) either the taxonomy tree, or the synth tree
 */

template<typename QT, typename TT, typename QM, typename TM>
json conflict_with_tree_impl(const QT & query_tree,
                             const TT & other_tree,
                             std::function<const QM*(const QM*,const QM*)> & query_mrca,
                             std::function<const TM*(const TM*,const TM*)> & other_mrca,
                             const witness_namer_t & witness_namer)
{
    conflict_stats stats(witness_namer);
    auto log_supported_by = [&stats](const QM* node2, const QM* node1) {
        if (is_fake_tip(node1)) {
            stats.add_terminal(node2,node1);
        } else {
            stats.add_supported_by(node2, node1);
        }
    };
    auto log_partial_path_of = [&stats](const QM* node2, const QM* node1) {
        if (is_fake_tip(node1)) {
            stats.add_terminal(node2,node1);
        } else {
            stats.add_partial_path_of(node2, node1);
        }
    };
    auto log_conflicts_with = [&stats](const QM* node2, const QM* node1) {
        stats.add_conflicts_with(node2, node1);
    };
    auto log_resolved_by = [&stats](const QM* node2, const QM* node1) {
        stats.add_resolved_by(node2, node1);
    };
    auto log_terminal = [&stats](const QM* node2, const QM* node1) {
        // Node1 might have an empty name if is a fake tip.
        if (not node1->get_name().empty()) {
            stats.add_terminal(node2, node1);
        }
    };

    {
        auto [induced_tree1, induced_tree2] = get_induced_trees2<ConflictTree>(query_tree, query_mrca, other_tree, other_mrca);

        perform_conflict_analysis(*induced_tree1,
                                  *induced_tree2,
                                  log_supported_by,
                                  log_partial_path_of,
                                  log_conflicts_with,
                                  log_resolved_by,
                                  log_terminal);
    }
/*
//  See PROBLEM notes above, on why we don't record when node2 resolves node1.
//    auto log_resolves = [&stats](const QM* node1, const QM* node2) {
//      if (not is_fake_tip(node1))
//          stats.add_resolves(node2, node1);
//    };
//    auto do_nothing = [](const QM*, const QM*) {};

    //The induced trees are destructively modified, so we can't reuse them.
    {
        auto induced_trees = get_induced_trees2<QT,TT,ConflictTree>(query_tree, query_mrca, other_tree, other_mrca);

        perform_conflict_analysis(*induced_trees.second,
                                  *induced_trees.first,
                                  do_nothing,
                                  do_nothing,
                                  do_nothing,
                                  log_resolves,
                                  do_nothing);
    }
*/

    {
        auto [induced_tree1, induced_tree2] = get_induced_trees2<ConflictTree>(query_tree, query_mrca, other_tree, other_mrca);
        return stats.get_json(*induced_tree1);
    }
}

json conflict_with_taxonomy(const ConflictTree& query_tree, const RichTaxonomy& Tax) {
    auto & taxonomy = Tax.get_tax_tree();

    using cfunc = std::function<const cnode_type*(const cnode_type*,const cnode_type*)>;
    using tfunc = std::function<const tnode_type*(const tnode_type*,const tnode_type*)>;
    cfunc query_mrca = [](const cnode_type* n1, const cnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };
    tfunc taxonomy_mrca = [](const tnode_type* n1, const tnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };
    witness_namer_t witness_namer = [&](const string& w) {return ott_witness_namer(w,Tax);};
    return conflict_with_tree_impl(query_tree, taxonomy, query_mrca, taxonomy_mrca, witness_namer);
}

json conflict_with_summary(const ConflictTree& query_tree,
                           const SummaryTree_t& summary,
                           const RichTaxonomy& Tax) {
    std::function<const cnode_type*(const cnode_type*,const cnode_type*)> query_mrca = [](const cnode_type* n1, const cnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };
    std::function<const snode_type*(const snode_type*,const snode_type*)> summary_mrca = [](const snode_type* n1, const snode_type* n2) {
        return mrca_from_depth(n1,n2);
    };
    witness_namer_t witness_namer = [&](const string& w) {return synth_witness_namer(w,summary,Tax);};
    return conflict_with_tree_impl(query_tree, summary, query_mrca, summary_mrca, witness_namer);
}

json conflict_with_newick(const ConflictTree& query_tree,
                          const ConflictTree& tree2,
                          const RichTaxonomy& Tax)
{
    std::function<const cnode_type*(const cnode_type*,const cnode_type*)> mrca = [](const cnode_type* n1, const cnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };
    witness_namer_t witness_namer = [&](const string& w) {return ott_witness_namer(w,Tax);};
    return conflict_with_tree_impl(query_tree, tree2, mrca, mrca, witness_namer);
}


unique_ptr<ConflictTree> get_induced_taxonomy(const ConflictTree& query_tree, const RichTaxonomy& taxonomy)
{
    using tfunc = std::function<const tnode_type*(const tnode_type*,const tnode_type*)>;

    tfunc taxonomy_mrca = [](const tnode_type* n1, const tnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };

    auto taxonomy_nodes_from_query_leaves = get_induced_nodes(query_tree, taxonomy.get_tax_tree());
    return get_induced_tree<ConflictTree>(taxonomy_nodes_from_query_leaves, taxonomy_mrca);
}

// Remove leaves (and their monotypic ancestors) in the query tree
// that are the ancestors of other leaves in the query tree.
void prune_ancestral_leaves(ConflictTree& query_tree, const RichTaxonomy& taxonomy)
{
    auto induced_taxonomy = get_induced_taxonomy(query_tree, taxonomy);

    auto ottid_to_induced_tax_node = get_ottid_to_node_map(*induced_taxonomy);

    LOG(WARNING)<<"induced taxonomy has "<<n_leaves(*induced_taxonomy)<<" leaves.";

    vector<cnode_type*> nodes_to_prune;
    for(auto leaf: iter_leaf(query_tree)) {
        auto ottid = leaf->get_ott_id();
        auto tax_node = ottid_to_induced_tax_node.at(ottid);
        if (not tax_node->is_tip()) {
            nodes_to_prune.push_back(leaf);
        }
    }
    for(auto node: nodes_to_prune) {
        delete_tip_and_monotypic_ancestors(query_tree, node);
    }
    LOG(WARNING)<<"query tree pruned down to "<<n_leaves(*induced_taxonomy)<<" leaves.";
}


// Remove leaves (and their monotypic ancestors) in the query tree
// that are the ancestors of other leaves in the query tree.
std::optional<OttId> has_ancestral_leaves(ConflictTree& query_tree, const RichTaxonomy& taxonomy)
{
    auto induced_taxonomy = get_induced_taxonomy(query_tree, taxonomy);

    auto ottid_to_induced_tax_node = get_ottid_to_node_map(*induced_taxonomy);

    vector<cnode_type*> nodes_to_prune;
    for(auto leaf: iter_leaf(query_tree)) {
        auto ottid = leaf->get_ott_id();
        auto tax_node = ottid_to_induced_tax_node.at(ottid);
        if (not tax_node->is_tip()) {
            return ottid;
        }
    }

    return {};
}

void check_and_forward_leaf_ott_ids(ConflictTree& query_tree, const RichTaxonomy& taxonomy, const string& tree_name)
{
    for(auto leaf: iter_leaf(query_tree))
    {
        if (not leaf->has_ott_id())
        {
            if (leaf->get_name().empty())
                throw OTCBadRequest()<<tree_name<<": Un-named leaf has no OTT id!";
            else
                throw OTCBadRequest()<<tree_name<<": Leaf '"<<leaf->get_name()<<"' has no OTT id!";
        }
        else
        {
            auto ottid = leaf->get_ott_id();
            auto valid_ottid = taxonomy.get_unforwarded_id(ottid);
            if (not valid_ottid)
            {
                if (leaf->get_name().empty())
                    throw OTCBadRequest()<<tree_name<<": Un-named leaf has bad OTT id "<<ottid<<"!";
                else
                    throw OTCBadRequest()<<tree_name<<": Leaf '"<<leaf->get_name()<<"' has bad OTT id "<<ottid<<"!";
            }
            else if (*valid_ottid != ottid)
                leaf->set_ott_id(*valid_ottid);
        }
    }
}

void check_all_nodes_have_node_names(const ConflictTree& query_tree, const string& tree_name) {
    for(const auto nd: iter_post_const(query_tree)) {
        if (nd->get_name().empty()) {
            auto E = OTCBadRequest();
            E << tree_name<<" has unnamed node";
            if (nd->has_ott_id()) {
                E << " with OTT Id=" << nd->get_ott_id();
            }
            E << "!";
            throw E;
        }
    }
}


// Find the smallest set C of taxa (leaf or internal) that we need to add as children of `id` so that
// i) all descendents of id are descendants of one of these children
// ii) all taxa in C are present in the summary tree
vector<OttId> extra_children_for_node(OttId id, const SummaryTree_t& summary, const RichTaxonomy& taxonomy) {
    // We could maybe improve speed by:
    // - following only nodes which are ancestors of synthesis leaves (i.e. unpruned leaves)
    // - following only nodes which are ancestors of input phylogeny leaves, exemplified to turn higher taxon leaves -> leaf leaves.
    auto& id_to_node = summary.get_data().id_to_node;
    auto& tax_id_to_node = taxonomy.get_tax_tree().get_data().id_to_node;
    // If the node is in synth already, then we don't need to add any children.
    if (id_to_node.count(id)) {
        return {};
    }
    // Growing list of descendants of `id` that are not in synth.
    vector<OttId> children;
    vector<OttId> bad_parents({id});
    // Walk a frontier leafward from id:
    for(std::size_t i=0;i<bad_parents.size();i++) {
        auto parent_id = bad_parents[i];
        auto tax_node = tax_id_to_node.at(bad_parents[i]);
        auto parent_node = taxonomy.included_taxon_from_id(parent_id);
        for(auto c: iter_child_const(*tax_node)) {
            auto child_id = c->get_ott_id();
            auto child_node = taxonomy.included_taxon_from_id(child_id);

            LOG(DEBUG) << "Lost taxon " << parent_node->get_name() << " ("<<parent_id<<") -> ";
            LOG(DEBUG) << "taxon " << child_node->get_name() << " (" << child_id << ") of degree " << child_node->get_out_degree() << ":";

            if (id_to_node.count(child_id)) {
                // In synth, we can stop expanding the frontier at this node.
                LOG(DEBUG) << "GOOD\n"; 
                children.push_back(child_id);
            } else {
                // Not in synth, we need to consider the children of this node.
                LOG(DEBUG) << "BAD\n";
                bad_parents.push_back(child_id);
            }
        }
    }
    return children;
}


// Find the smallest set C of taxa (leaf or internal) that we need to add as children of `id` so that
// i) all descendents of id are descendants of one of these children
// ii) all taxa in C are present in the summary tree
map<cnode_type*, vector<OttId>>
extra_children_for_tree(ConflictTree& query_tree, const ConflictTree& other_tree, const RichTaxonomy& taxonomy)
{
    auto induced_taxonomy = get_induced_taxonomy(other_tree, taxonomy);

    auto ottid_to_induced_tax_node = get_ottid_to_const_node_map(*induced_taxonomy);

    // For each leaf of the query tree
    map<cnode_type*, vector<OttId>> leaves_to_add;
    for(auto leaf: iter_leaf(query_tree))
    {
        auto leaf_id = leaf->get_ott_id();

        // .. that is in the induced taxonomy ...
        if (not ottid_to_induced_tax_node.count(leaf_id)) continue;

        // .. but is not a tip ...
        auto tax_node = ottid_to_induced_tax_node.at(leaf_id);

        if (tax_node->is_tip()) continue;

        // ... find all the descendant tips in the other tree.
        vector<OttId> children;
        for(auto leaf2: iter_post_n(*tax_node))
            if (leaf2->is_tip())
                children.push_back(leaf2->get_ott_id());

        leaves_to_add.insert({leaf, children});
    }

    return leaves_to_add;
}

void add_children(ConflictTree& tree, const map<cnode_type*, vector<OttId>>& children_to_add)
{
    for(const auto& [leaf, child_ids]: children_to_add)
        for(auto id: child_ids)
        {
            auto new_leaf = tree.create_child(leaf);
            new_leaf->set_ott_id(id);
            // All nodes must have a name!
            new_leaf->set_name("ott"+std::to_string(id));
        }
}

/*
 * OK, a general approach to normalizing two trees for comparison:
 * 0. All tips must have OTT IDs.
 *    The taxonomy and ids are used to assess if a tip is identical, ancestral, or descended from another tip.
 * 1. Remove duplicate tips.
 * 2. Remove tips that are ancestral to other tips.
 * At this point, we COULD just add as children all OTT leaves that are descendants of each higher taxon.
 *
 * A possibly simpler approach would be to.
 * 3. Remove tips that have no shared descendants with tips in the other tree.
 * 4. For each tip, add as children tips from the other tree that are descendants.
 */

string conflict_ws_method(const SummaryTree_t& summary,
                          const RichTaxonomy & taxonomy,
                          std::unique_ptr<ConflictTree>& query_tree,
                          const string& tree2s)
{
    // 0. Prune unmapped leaves -- this already handles forwards.
    auto leaf_counts = prune_unmapped_leaves(*query_tree, taxonomy);
    if (leaf_counts.first < 3) {
        throw OTCBadRequest()<<"Query tree has only "<<leaf_counts.first<<" leaves with an OTT id!";
    }
    // 1. Check that all leaves in input tree have OTT ids
    check_and_forward_leaf_ott_ids(*query_tree, taxonomy, "tree1");
    // 2. Check that all leaves in input tree have node names
    check_all_nodes_have_node_names(*query_tree, "tree1");
    // 3. Prune leaves with duplicate ott ids
    prune_duplicate_ottids(*query_tree);
    // 4. Prune leaves of the query that are ancestral to other query leaves
    prune_ancestral_leaves(*query_tree, taxonomy);
    // 5. Add children to higher-taxon leaves of query_tree that are not in the synth tree.
    if (tree2s == "synth") {
        map<cnode_type*,vector<OttId>> children_to_add;
        for(auto leaf: iter_leaf(*query_tree)) {
//          LOG(WARNING)<<"Considering leaf "<<leaf->get_ott_id()<<":";
            auto leaf_id = leaf->get_ott_id();
            auto c = extra_children_for_node(leaf_id, summary, taxonomy);
            if (not c.empty()) {
                children_to_add.insert({leaf,c});
                LOG(DEBUG)<<"ottid "<<leaf->get_ott_id()<<" has frontier ";
                for(auto id: c) {
                    LOG(DEBUG)<<"   "<<id;
                }
            }
        }
        // Add nodes with the specified ottids
        add_children(*query_tree, children_to_add);
    }
    if (tree2s == "ott")
    {
        // Compute depth and number of tips for query_tree.
        compute_depth(*query_tree);
        compute_tips(*query_tree);

        return conflict_with_taxonomy(*query_tree, taxonomy).dump(1);
    }
    else if (tree2s == "synth")
    {
        // Compute depth and number of tips for query_tree.
        compute_depth(*query_tree);
        compute_tips(*query_tree);

        return conflict_with_summary(*query_tree, summary, taxonomy).dump(1);
    }
    else if (tree2s.size() > 0 and tree2s[0] == '(') {
        auto tree2 = tree_from_newick_string<ConflictTree>(tree2s);

        // 0. Prune unmapped leaves -- this already handles forwards.
        auto leaf_counts = prune_unmapped_leaves(*tree2, taxonomy);
        if (leaf_counts.first < 3) {
            throw OTCBadRequest()<<"tree2 tree has only "<<leaf_counts.first<<" leaves with an OTT id!";
        }
        // 1. Check that all leaves in input tree have OTT ids
        check_and_forward_leaf_ott_ids(*tree2, taxonomy, "tree2");
        // 2. Check that all leaves in tree2 have node names
        check_all_nodes_have_node_names(*tree2, "tree2");
        // 3. Prune leaves with duplicate ott ids
        prune_duplicate_ottids(*tree2);
        // 4. Prune leaves of the query that are ancestral to other query leaves
        prune_ancestral_leaves(*tree2, taxonomy);

        // Find children in each tree that are descendents of nodes in the other tree.
        // Get the children before modifying either tree.
        auto children_to_add1 = extra_children_for_tree(*query_tree, *tree2, taxonomy);
        auto children_to_add2 = extra_children_for_tree(*tree2, *query_tree, taxonomy);

        // Add nodes with the specified ottids
        add_children(*query_tree, children_to_add1);
        add_children(*tree2, children_to_add2);

        // Compute depth and number of tips for query_tree.
        compute_depth(*query_tree);
        compute_tips(*query_tree);

        // Compute depth and number of tips for tree2
        compute_depth(*tree2);
        compute_tips(*tree2);

        return conflict_with_newick(*query_tree, *tree2, taxonomy).dump(1);
    }
    throw OTCBadRequest() << "tree2 = '" << tree2s << "' not recognized!";
}

string newick_conflict_ws_method(const SummaryTree_t& summary,
                                 const RichTaxonomy & taxonomy,
                                 const string& tree1s,
                                 const string& tree2s) {
    try {
        LOG(WARNING)<<"newick conflict:";
        LOG(DEBUG)  <<"  tree1s = "<<tree1s;
        LOG(DEBUG)  <<"  tree2s = "<<tree2s;
        auto query_tree = tree_from_newick_string<ConflictTree>(tree1s);
        return conflict_ws_method(summary, taxonomy, query_tree, tree2s);
    } catch (otc::OTCParsingError& e) {
        LOG(WARNING)<<"newick conflict: error parsing newick string '"<<tree1s<<"'";
        throw OTCBadRequest() << "Error parsing newick:  \n:"<<e.what();
    }
}

string phylesystem_conflict_ws_method(const SummaryTree_t& summary,
                                      const RichTaxonomy & taxonomy,
                                      const string& tree1s,
                                      const string& tree2s) {
    LOG(WARNING)<<"phylesystem conflict:";
    LOG(DEBUG)  <<"  tree1s = "<<tree1s;
    LOG(DEBUG)<<"  tree2s = "<<tree2s;
    auto query_tree = get_phylesystem_tree<ConflictTree>(tree1s);
    return conflict_ws_method(summary, taxonomy, query_tree, tree2s);
}


} // namespace otc
