#include <regex>
#include "otc/ws/tolws.h"
#include "otc/ws/tolwsadaptors.h"
#include "otc/ws/trees_to_serve.h"
#include "otc/conflict.h"
#include "otc/ws/nexson/nexson.h"
#include "otc/ws/prune.h"
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



// Get the subtree of T1 connecting the leaves of T1 that are also in T2.
template <typename Tree_Out_t, typename Tree1, typename Tree2>
pair<unique_ptr<Tree_Out_t>,unique_ptr<Tree_Out_t>>
get_induced_trees2(Tree1& T1,
                   std::function<const typename Tree1::node_type*(const typename Tree1::node_type*,const typename Tree1::node_type*)> MRCA_of_pair1,
                   Tree2& T2,
                   std::function<const typename Tree2::node_type*(const typename Tree2::node_type*,const typename Tree2::node_type*)> MRCA_of_pair2)
{
    LOG(WARNING)<<"T1 = "<<newick_string(T1);
    LOG(WARNING)<<"n_leaves(T1) = "<<n_leaves(T1);
    LOG(WARNING)<<"n_leaves(T2) = "<<n_leaves(T2);

    // 1. First construct the induced tree for T2.

    // 1a. Find the nodes of T2 that corresponds to leaves of T1.
    //     Note that some of these nodes could be ancestral to other ones in T2.
    auto T2_nodes_from_T1_leaves = get_induced_nodes(T1,T2);
    LOG(WARNING)<<T2_nodes_from_T1_leaves.size()<<" leaves of T1 are in T2";

    // 1b. Actually construct the induced tree for T2.
    //     It might have fewer leaves than in T2_nodes_from_T1_leaves, if some of the nodes are ancestral to others.
    auto induced_tree2 = get_induced_tree<Tree2, Tree_Out_t>(T2_nodes_from_T1_leaves, MRCA_of_pair2);
    LOG(WARNING)<<"n_leaves(induced_tree2) = "<<n_leaves(*induced_tree2);
    LOG(WARNING)<<"induced-tree2a = "<<newick_string(*induced_tree2);

    // 1c. Rename internal nodes of synth/taxonomy to ottXXX instead of taxon name
    for(auto node: iter_post(*induced_tree2)) {
        if (node->has_ott_id()) {
            node->set_name("ott"+std::to_string(node->get_ott_id()));
        }
    }
    LOG(WARNING)<<"induced_tree2 = "<<newick_string(*induced_tree2);

    //FIXME - handle cases like (Homo sapiens, Homo) by deleting monotypic nodes at the root.

    // 2. Keep only leaves from T1 that (a) have an ottid in T2 and (b) map to a leaf of induced_tree2
    auto ottid_to_induced_tree2_node = otc::get_ottid_to_const_node_map(*induced_tree2);
    std::vector<const typename Tree1::node_type*> induced_leaves1;
    for(auto leaf: iter_leaf_const(T1))
    {
        if (not leaf->has_ott_id())
        {
            LOG(WARNING)<<"Dropping tip: no ott id";
            continue;
        }

        auto it = ottid_to_induced_tree2_node.find(leaf->get_ott_id());
        if (it == ottid_to_induced_tree2_node.end()) {
            LOG(WARNING)<<"Dropping tip "<<leaf->get_ott_id()<<": not found in induced taxonomy-or-synth tree.";
        } else if (not it->second->is_tip()) {
            LOG(WARNING)<<"Dropping higher taxon tip "<<leaf->get_ott_id();
        } else {
            induced_leaves1.push_back(leaf);
        }
    }
    LOG(WARNING)<<"keeping "<<induced_leaves1.size()<<" leaves from T1";

    // 3. Construct the induced tree for T1
    auto induced_tree1 = get_induced_tree<Tree1, Tree_Out_t>(induced_leaves1, MRCA_of_pair1);
    LOG(WARNING) << "n_leaves(induced_tree1) = " << n_leaves(*induced_tree1);
    LOG(WARNING) << "induced_tree1 = " << newick_string(*induced_tree1);
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

struct conflict_stats {
    protected:
    template <typename T>
    void throw_otcerror_if_second_not_true(const T & result, const cnode_type* node1) {
        if (not result.second) {
            throw OTCError() << "key "<<node_name(node1) << " occurs twice in a conflict_stats map!";
        }    
    }
    public:
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
    json get_json(const ConflictTree&, const RichTaxonomy&) const;
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

inline json get_node_status(const string& witness, string status, const RichTaxonomy& Tax) {
    json j;
    j["witness"] = extract_node_name_if_present(witness);
    j["status"] = std::move(status);
    long raw_ott_id = long_ott_id_from_name(witness);
    if (raw_ott_id >= 0) {
        OttId id = check_ott_id_size(raw_ott_id);
        auto nd = Tax.included_taxon_from_id(id);
        if (nd != nullptr) {
            j["witness_name"] = nd->get_name();
        }
    }
    return j;
}

json get_conflict_node_status(const set<pair<string,int>>& witnesses, string status, const RichTaxonomy& Tax)
{
    json j_witnesses = json::array();
    json j_witness_names = json::array();
    for(auto& [node_name, _]: witnesses)
    {
        json witness_name;

        if (long raw_ott_id = long_ott_id_from_name(node_name); raw_ott_id >= 0)
        {
            OttId id = check_ott_id_size(raw_ott_id);
            if (auto nd = Tax.included_taxon_from_id(id))
                witness_name  = nd->get_name();
        }
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
// Let's solve this situation by NOT reporting when tree2 resolves tree1.


json conflict_stats::get_json(const ConflictTree& tree, const RichTaxonomy& Tax) const {
    json nodes;
//    for(auto& x: resolves) {
//        nodes[extract_node_name_if_present(x.first)] = get_node_status(x.second, "resolves", Tax);
//    }
    for(auto& x: resolved_by) {
        nodes[extract_node_name_if_present(x.first)] = get_node_status(x.second, "resolved_by", Tax);
    }
    for(auto& x: supported_by) {
        nodes[extract_node_name_if_present(x.first)] = get_node_status(x.second, "supported_by", Tax);
    }
    for(auto& x: partial_path_of) {
        nodes[extract_node_name_if_present(x.first)] = get_node_status(x.second, "partial_path_of", Tax);
    }
    for(auto& x: terminal) {
        nodes[extract_node_name_if_present(x.first)] = get_node_status(x.second, "terminal", Tax);
    }
    for(auto& x: conflicts_with) {
        nodes[extract_node_name_if_present(x.first)] = get_conflict_node_status(x.second, "conflicts_with", Tax);
    }
    // For monotypic nodes in the query, copy annotation from child.
    for(auto it: iter_post_const(tree)) {
        if (it->is_outdegree_one_node()) {
            auto name = extract_node_name_if_present(it->get_name());
            auto child_name = extract_node_name_if_present(it->get_first_child()->get_name());
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
                             const RichTaxonomy& Tax)
{
    conflict_stats stats;
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
        auto induced_trees = get_induced_trees2<ConflictTree>(query_tree, query_mrca, other_tree, other_mrca);

        perform_conflict_analysis(*induced_trees.first,
                                  *induced_trees.second,
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
        auto induced_trees = get_induced_trees2<ConflictTree>(query_tree, query_mrca, other_tree, other_mrca);
        return stats.get_json(*induced_trees.first, Tax);
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
    return conflict_with_tree_impl(query_tree, taxonomy, query_mrca, taxonomy_mrca, Tax);
}

json conflict_with_summary(const ConflictTree& query_tree,
                           const SummaryTree_t& summary,
                           const RichTaxonomy& Tax) {
    using snode_type = SummaryTree_t::node_type;
    std::function<const cnode_type*(const cnode_type*,const cnode_type*)> query_mrca = [](const cnode_type* n1, const cnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };
    std::function<const snode_type*(const snode_type*,const snode_type*)> summary_mrca = [](const snode_type* n1, const snode_type* n2) {
        return find_mrca_via_traversal_indices(n1,n2);
    };
    return conflict_with_tree_impl(query_tree, summary, query_mrca, summary_mrca, Tax);
}


// Remove leaves (and their monotypic ancestors) in the query tree
// that are the ancestors of other leaves in the query tree.
void prune_ancestral_leaves(ConflictTree& query_tree, const RichTaxonomy& taxonomy) {
    using tfunc = std::function<const tnode_type*(const tnode_type*,const tnode_type*)>;
    tfunc taxonomy_mrca = [](const tnode_type* n1, const tnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };

    auto taxonomy_nodes_from_query_leaves = get_induced_nodes(query_tree, taxonomy.get_tax_tree());
    auto induced_taxonomy = get_induced_tree<const RichTaxTree, ConflictTree>(taxonomy_nodes_from_query_leaves,
                                                                              taxonomy_mrca);
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

void check_all_leaves_have_ott_ids(const ConflictTree& query_tree) {
    for(auto leaf: iter_leaf_const(query_tree)) {
        if (leaf->has_ott_id()) {
            continue;
        }
        if (leaf->get_name().empty()) {
            throw OTCBadRequest()<<"Un-named leaf has no OTT id!";
        } else {
            throw OTCBadRequest()<<"Leaf '"<<leaf->get_name()<<"' has no OTT id!";
        }
    }
}

void check_all_nodes_have_node_names(const ConflictTree& query_tree) {
    for(const auto nd: iter_post_const(query_tree)) {
        if (nd->get_name().empty()) {
            auto E = OTCBadRequest();
            E << "Query tree has unnamed node";
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
    // 0. Prune unmapped leaves.
    auto leaf_counts = prune_unmapped_leaves(*query_tree, taxonomy);
    if (leaf_counts.first < 3) {
        throw OTCBadRequest()<<"Query tree has only "<<leaf_counts.first<<" leaves with an OTT id!";
    }
    // 1. Check that all leaves in input tree have OTT ids
    check_all_leaves_have_ott_ids(*query_tree);
    // 2. Check that all leaves in input tree have node names
    check_all_nodes_have_node_names(*query_tree);
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
                LOG(WARNING)<<"ottid "<<leaf->get_ott_id()<<" has frontier ";
                for(auto id: c) {
                    LOG(WARNING)<<"   "<<id;
                }
            }
        }
        // Add nodes with the specified ottids
        for(const auto& job: children_to_add) {
            auto leaf = job.first;
            auto& child_ids = job.second;
            for(auto id: child_ids) {
                query_tree->create_child(leaf)->set_ott_id(id);
            }
        }
    }
    // 5. Compute depth and number of tips for query_tree.
    compute_depth(*query_tree);
    compute_tips(*query_tree);
    if (tree2s == "ott") {
        return conflict_with_taxonomy(*query_tree, taxonomy).dump(1);
    } else if (tree2s == "synth") {
        return conflict_with_summary(*query_tree, summary, taxonomy).dump(1);
    }
    throw OTCBadRequest() << "tree2 = '" << tree2s << "' not recognized!";
}

string newick_conflict_ws_method(const SummaryTree_t& summary,
                                 const RichTaxonomy & taxonomy,
                                 const string& tree1s,
                                 const string& tree2s) {
    try {
        LOG(WARNING)<<"newick conflict: tree1s = '"<<tree1s<<"'   tree1s = '"<<tree2s<<"'";
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
    LOG(WARNING)<<"phylesystem conflict: tree1s = '"<<tree1s<<"'   tree1s = '"<<tree2s<<"'";
    auto query_tree = get_phylesystem_tree<ConflictTree>(tree1s);
    return conflict_ws_method(summary, taxonomy, query_tree, tree2s);
}


} // namespace otc
