// See https://github.com/OpenTreeOfLife/opentree/wiki/Open-Tree-of-Life-APIs-v3#synthetic-tree
#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include <tuple>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include "json.hpp"

using namespace otc;
using json = nlohmann::json;
using std::vector;
template <typename T, typename U>
using Map = std::unordered_multimap<T,U>;
template <typename T, typename U>
using map = std::unordered_map<T,U>;
template <typename T>
using set = std::unordered_set<T>;
using std::string;
//using std::map;
using std::pair;
using std::tuple;

// TODO: Could we exemplify tips here, if we had access to the taxonomy?

namespace std
{
template<typename T, typename U>
struct hash<pair<T,U>>
{
    std::size_t operator()(const std::pair<T,U>& p) const noexcept {return std::hash<T>()(p.first) * std::hash<U>()(p.second);}
};
}

struct RTNodeDepth {
    int depth = 0; // depth = number of nodes to the root of the tree including the  endpoints (so depth of root = 1)
    int n_tips = 0;
    int n_include_tips = 0;
    RootedTreeNode<RTNodeDepth>* summary_node;
};

using Tree_t = RootedTree<RTNodeDepth, RTreeNoData>;
using node_t = Tree_t::node_type;

int n_include_tips(const node_t* node)
{
    return node->getData().n_include_tips;
}

int n_tips(const node_t* node)
{
    return node->getData().n_tips;
}

int& n_include_tips(node_t* node)
{
    return node->getData().n_include_tips;
}

int& n_tips(node_t* node)
{
    return node->getData().n_tips;
}

int depth(const Tree_t::node_type* node);
int& depth(Tree_t::node_type* node);
Tree_t::node_type* summary_node(const Tree_t::node_type* node);
Tree_t::node_type*& summary_node(Tree_t::node_type* node);
Tree_t::node_type* nmParent(Tree_t::node_type* node);
void computeSummaryLeaves(Tree_t& tree, const map<long,Tree_t::node_type*>& summaryOttIdToNode);
string getSourceNodeNameIfAvailable(const Tree_t::node_type* node);
Tree_t::node_type* get_root(Tree_t::node_type* node);
const Tree_t::node_type* get_root(const Tree_t::node_type* node);
void find_anc_conflicts(Tree_t::node_type* node, vector<Tree_t::node_type*>& conflicts);
void find_conflicts(const Tree_t& tree, vector<Tree_t::node_type*>& conflicts);

inline int depth(const Tree_t::node_type* node) {
    assert(node->getData().depth > 0);
    return node->getData().depth;
}


inline int& depth(Tree_t::node_type* node) {
    assert(node->getData().depth > 0);
    return node->getData().depth;
}

inline Tree_t::node_type* summary_node(const Tree_t::node_type* node) {
    return node->getData().summary_node;
}

inline Tree_t::node_type*& summary_node(Tree_t::node_type* node) {
    return node->getData().summary_node;
}

// returns the most recent non-monotypic ancestor of node (or nullptr if there is no such node)
inline Tree_t::node_type* nmParent(Tree_t::node_type* node) {
    do {
        node = node->getParent();
    } while (node and node->isOutDegreeOneNode());
    return node;
}

// Compute number of tips at this node or below it
void compute_tips(Tree_t& tree)
{
    // Iterate over nodes, leaves first
    for(auto nd: iter_post(tree))
    {
        if (nd->isTip())
            nd->getData().n_tips = 1;

        auto p = nd->getParent();
        if (p)
            p->getData().n_tips += nd->getData().n_tips;
    }
}

// uses the OTT Ids in `tree` to fill in the `summary_node` field of each leaf
void computeSummaryLeaves(Tree_t& tree, const map<long,Tree_t::node_type*>& summaryOttIdToNode) {
    for(auto leaf: iter_leaf(tree)) {
        summary_node(leaf) = summaryOttIdToNode.at(leaf->getOttId());
    }
}


string getSourceNodeNameIfAvailable(const Tree_t::node_type* node) {
    string name = node->getName();
    auto source = getSourceNodeName(name);
    if (source)
        return *source;
    else
        return name;
}

Tree_t::node_type* get_root(Tree_t::node_type* node) {
    while (node->getParent()) {
        node = node->getParent();
    }
    return node;
}

const Tree_t::node_type* get_root(const Tree_t::node_type* node)
{
    while (node->getParent()){
        node = node->getParent();
    }
    return node;
}

// assumes `node` is marked with `bits`. Returns the parent of `node` and
//  also marks the parent with `bits` as a side effect.
const node_t* trace_to_parent(const node_t* node, set<const node_t*>& nodes) {
    assert(nodes.count(node));
    node = node->getParent(); // move to parent and insert it.
    nodes.insert(node);
    return node;
}

node_t* trace_to_parent(node_t* node, set<node_t*>& nodes) {
    assert(nodes.count(node));
    node = node->getParent(); // move to parent and insert it.
    nodes.insert(node);
    return node;
}

/// Walk up the tree from node1 and node2 until we find the common ancestor, putting nodes into set `nodes`.
template <typename N>
N* trace_find_MRCA(N* node1, N* node2, set<N*>& nodes)
{
    assert(node1 or node2);
    if (not node1) {
        assert(nodes.count(node2));
        return node2;
    }
    if (not node2) {
        assert(nodes.count(node2));
        return node1;
    }
    assert(node1 and node2);
    assert(get_root(node1) == get_root(node2));
    assert(nodes.count(node1));
    assert(nodes.count(node2));
    while (depth(node1) > depth(node2)){
        node1 = trace_to_parent(node1, nodes);
    }
    while (depth(node1) < depth(node2)){
        node2 = trace_to_parent(node2, nodes);
    }
    assert(depth(node1) == depth(node2));
    while (node1 != node2) {
        assert(node1->getParent());
        assert(node2->getParent());
        node1 = trace_to_parent(node1, nodes);
        node2 = trace_to_parent(node2, nodes);
    }
    assert(node1 == node2);
    assert(nodes.count(node1));
    assert(nodes.count(node2));
    return node1;
}

template <typename N>
vector<N*> leaf_nodes_below(N* node)
{
    vector<N*> nodes;
    for(auto nd: iter_leaf_n_const(*node))
        nodes.push_back(nd);
    return nodes;
}

vector<node_t*> map_to_summary(const vector<const node_t*>& nodes)
{
    vector<node_t*> nodes2(nodes.size(),nullptr);
    for(int i=0;i<nodes2.size();i++)
        nodes2[i] = summary_node(nodes[i]);
    return nodes2;
}

auto find_induced_nodes(const vector<node_t*>& leaves)
{
    set<node_t*> nodes;
    node_t* MRCA = nullptr;
    for(auto leaf: leaves)
    {
        nodes.insert(leaf);
        MRCA = trace_find_MRCA(MRCA, leaf, nodes);
    }
    vector<node_t*> vnodes;
    for(auto nd: nodes)
        vnodes.push_back(nd);
    return vnodes;
}

std::unique_ptr<Tree_t> get_induced_tree(const std::vector<const node_t*>& leaves)
{
    // 1. Find MRCA and all nodes in the tree
    const Tree_t::node_type* MRCA = nullptr;
    set<const node_t*> nodes;
    for(auto leaf: leaves)
    {
        nodes.insert(leaf);
        MRCA = trace_find_MRCA(MRCA, leaf, nodes);
    }

    std::unique_ptr<Tree_t> induced_tree(new Tree_t());

    // 2. Construct duplicate nodes for the induced tree, recording correspondence
    map<const node_t*, node_t*> to_induced_tree;
    for(auto nd: nodes)
    {
        auto nd2 = induced_tree->createNode(nullptr);

        if (nd->hasOttId())
            nd2->setOttId(nd->getOttId());

        if (nd->getName().size())
            nd2->setName(nd->getName());

        to_induced_tree[nd] = nd2;
    }

    // 3. Link corresponding nodes to their corresponding parents
    for(auto nd: nodes)
    {
        auto p = nd->getParent();

        auto nd2 = to_induced_tree.find(nd)->second;
        auto p2_it = to_induced_tree.find(p);

        if (p2_it == to_induced_tree.end())
        {
            assert(nd == MRCA);
        }
        else
        {
            auto p2 = p2_it->second;
            p2->addChild(nd2);
        }
    }

    // 4. Set the root of the induced tree to node corresponding to the MRCA
    induced_tree->_setRoot( to_induced_tree.at(MRCA) );
    
    return induced_tree;
}

// Get a list of leaves of tree 1 that are also in tree 2.
vector<const Tree_t::node_type*> get_induced_leaves(const Tree_t& T1, const map<long, const Tree_t::node_type*>& nodes1,
                                                    const Tree_t& T2, const map<long, const Tree_t::node_type*>& nodes2)
{
    vector<const Tree_t::node_type*> leaves;
    
    if (nodes2.size() < nodes1.size())
    {
        for(auto leaf: iter_leaf_const(T2))
        {
            auto id = leaf->getOttId();
            auto it = nodes1.find(id);
            if (it != nodes1.end())
                leaves.push_back(it->second);
        }
    }
    else
    {
        for(auto leaf: iter_leaf_const(T1))
        {
            auto id = leaf->getOttId();
            if (nodes2.find(id) != nodes2.end())
                leaves.push_back(leaf);
        }
    }

    return leaves;
}

// Get the subtree of T1 connecting the leaves of T1 that are also in T2.
std::unique_ptr<Tree_t> get_induced_tree(const Tree_t& T1, const map<long, const Tree_t::node_type*>& nodes1,
                                         const Tree_t& T2, const map<long, const Tree_t::node_type*>& nodes2)
{
    auto induced_leaves = get_induced_leaves(T1, nodes1, T2, nodes2);

    auto induced_tree =  get_induced_tree(induced_leaves);
    induced_tree->setName(T1.getName());

    return induced_tree;
}


json get_support_blob_as_array(const Map<string,string>& M)
{
    json support_blob = json::object();
    auto m_it = M.begin();
    for (;  m_it != M.end();  )
    {
        string source_id = (*m_it).first;
        
        json nodes_array = json::array();
        
        // Iterate over all map elements with key == source_id
        auto keyRange = M.equal_range(source_id);
        for (auto s_it = keyRange.first;  s_it != keyRange.second;  ++s_it)
            nodes_array.push_back((*s_it).second);
        
        support_blob[source_id] = nodes_array;

        m_it = keyRange.second;
    }
    return support_blob;
}

json get_support_blob_as_single_element(const Map<string,string>& M)
{
    json support_blob = json::object();
    auto m_it = M.begin();
    for (;  m_it != M.end();  )
    {
        string source_id = (*m_it).first;
        
        json element;
        
        // Iterate over all map elements with key == source_id
        auto keyRange = M.equal_range(source_id);
        int count = 0;
        for (auto s_it = keyRange.first;  s_it != keyRange.second;  ++s_it)
        {
            element = (*s_it).second;
            ++count;
        }
        assert(count == 1);

        support_blob[source_id] = element;

        m_it = keyRange.second;
    }
    return support_blob;
}

void set_support_blob_as_array(json& j, const map<string,Map<string,string>>& m, const string& field, const string& name)
{
    json support_blob;
    auto it = m.find(name);
    if (it != m.end())
        support_blob = get_support_blob_as_array(it->second);
    if (not support_blob.empty())
        j[field] = support_blob;
}

void set_support_blob_as_single_element(json& j, const map<string,Map<string,string>>& m, const string& field, const string& name)
{
    json support_blob;
    auto it = m.find(name);
    if (it != m.end())
        support_blob = get_support_blob_as_single_element(it->second);
    if (not support_blob.empty())
        j[field] = support_blob;
}

void add_element(map<string, Map<string, string>>& m, map<string, set<pair<string,string>>>& s,
                 const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const string& source)
{
    string synth = synth_node->getName();
    string node = getSourceNodeNameIfAvailable(input_node);

    pair<string,string> x{source, node};

    if (not s[synth].count(x))
    {
        s[synth].insert(x);
        m[synth].insert(x);
    }
}

map<long, const node_t*> get_ottid_to_const_node_map(const Tree_t& T)
{
    map<long, const node_t*> ottid;
    for(auto nd: iter_pre_const(T))
        if (nd->hasOttId())
            ottid[nd->getOttId()] = nd;

    return ottid;
}

map<long, node_t*> get_ottid_to_node_map(Tree_t& T)
{
    map<long, node_t*> ottid;
    for(auto nd: iter_pre(T))
        if (nd->hasOttId())
            ottid[nd->getOttId()] = nd;

    return ottid;
}

void remove_monotypic_node(node_t* nd)
{
    auto child = nd->getFirstChild();
    child->detachThisNode();
    nd->addSibOnRight(child);
    nd->detachThisNode();
}

map<string,string> suppressAndRecordMonotypic(Tree_t& tree)
{
    map<string,string> to_child;

    std::vector<Tree_t::node_type*> remove;
    for (auto nd:iter_post(tree)) {
        if (nd->isOutDegreeOneNode()) {
            remove.push_back(nd);
        }
    }

    for (auto nd: remove)
    {
        assert(nd->isOutDegreeOneNode());
        auto child = nd->getFirstChild();
        assert(not child->isOutDegreeOneNode());

        if (nd->getName().size())
        {
            assert(to_child.count(nd->getName()) == 0);
            to_child[nd->getName()] = child->getName();
            remove_monotypic_node(nd);
            delete nd;
        }
    }
    return to_child;
}

void destroy_children(node_t* node)
{
    vector<node_t*> nodes;
    while(auto n = node->getFirstChild())
    {
        n->detachThisNode();
        nodes.push_back(n);
    }

    for(std::size_t i = 0; i < nodes.size(); i++) {
        while(auto n = node->getFirstChild())
        {
            n->detachThisNode();
            nodes.push_back(n);
        }
        delete nodes[i];
    }

    assert(node->isTip());
}

typedef std::function<void(const node_t* node2, const node_t* node1)> node_logger_t;

// Map nodes of T1 onto T2.
// * Monotypic nodes of T1 are ignored.
// * Leaves of T1 should correspond to leaves of T2, or nothing in T2. (?)
// * Each non-monotypic node of T1 will equal (supported_by/partial_path_of), conflict with, or be a terminal of nodes in T2.
// * A node of T1 can equal one node of T2 (supported_by) or several nodes of T2 (partial_path_of).
// * A node of T1 can conflict with several nodes of T2.
// * A node of T1
void perform_conflict_analysis(const Tree_t& tree1,
                               const map<long, const node_t*>& ottid_to_node1,
                               const Tree_t& tree2,
                               const map<long, const node_t*>& ottid_to_node2,
                               node_logger_t log_supported_by,
                               node_logger_t log_partial_path_of,
                               node_logger_t log_conflicts_with,
                               node_logger_t log_resolved_by,
                               node_logger_t log_terminal)
{
    // Handle non-leaf correspondence?
    auto induced_tree1 = get_induced_tree(tree1, ottid_to_node1, tree2, ottid_to_node2);
    auto induced_tree2 = get_induced_tree(tree2, ottid_to_node2, tree1, ottid_to_node1);

    computeDepth(*induced_tree1);
    compute_tips(*induced_tree1);

    computeDepth(*induced_tree2);
    compute_tips(*induced_tree2);

    vector<Tree_t::node_type*> conflicts;

    auto map1 = get_ottid_to_node_map(*induced_tree1);
    auto map2 = get_ottid_to_node_map(*induced_tree2);
        
    // make summary_node field of induced_tree1 leaves point to leaves of induced_tree2.
    for(auto leaf: iter_leaf(*induced_tree1))
    {
        auto leaf2 = map2.at(leaf->getOttId());
        summary_node(leaf) = leaf2;
    }
        
    int L = countLeaves(*induced_tree1);
    assert(L == countLeaves(*induced_tree2));
        
    vector<node_t*> tree_nodes;
    for(auto nd: iter_post(*induced_tree1))
        tree_nodes.push_back(nd);

    for(auto nd: tree_nodes)
    {
        if (not nd->getParent()) continue;

        // Ignore knuckles in input trees.
        // (Note that in general, if we've pruned this tree down to match the shared taxon set
        //  then this could produce knuckles that were not originally there.)
        if (nd->isOutDegreeOneNode()) continue;

        // If this node contains all tips under it, then it doesn't correspond to a split.
        if (nd->getData().n_tips == L) continue;
            
        // If this node is a tip, the mark the corresponding nodes
        if (nd->isTip())
        {
            auto nd2 = summary_node(nd);
            log_terminal(nd2, nd);
            nd2 = nd2->getParent();
            for(;nd2 and nd2->isOutDegreeOneNode();nd2 = nd2->getParent())
                log_terminal(nd2, nd);
            continue;
        }

        // Find the list of nodes in the input tree that are below nd.
        auto leaves1 = leaf_nodes_below(const_cast<const node_t*>(nd));

        // Since nd is not a tip, and not monotypic, it should have at least 2 leaves below it.
        assert(leaves1.size() >= 2);

        int L2 = 0;
        for(auto nd: leaves1)
            L2 += n_tips(nd);

        // Find the corresponding list of nodes in the summary tree
        auto leaves2 = map_to_summary(leaves1);

        // Find the nodes in the induced tree of those nodes
        vector<node_t*> nodes = find_induced_nodes(leaves2);

        // Sort the nodes by depth to ensure all children occur before their parents -- unfortunately n*log(n) !
        std::sort(nodes.begin(), nodes.end(), [](auto x, auto y){return depth(x) > depth(y);});

        // The MRCA should be the last node in the vector.
        auto MRCA = nodes.back();

        // The n_include_tips for a parent node should count the n_include_tips for this node
        for(int i=0;i<nodes.size()-1;i++)
        {
            auto nd = nodes[i];
            if (nd->isTip())
                n_include_tips(nd) = n_tips(nd);
            auto p = nd->getParent();
            assert(p);
            assert(nd != MRCA);
            n_include_tips(p) += n_include_tips(nd);
            assert(n_include_tips(nd) <= n_tips(nd));
        }
            
        // If MRCA includes all and only the tips under nd, then MRCA is supporting or partial_path_of
        bool conflicts_or_resolved_by = n_include_tips(MRCA) < n_tips(MRCA);

        // Supported_by or partial_path_of
        if (not conflicts_or_resolved_by)
        {
            assert(MRCA->getParent());
            if (MRCA->getParent()->getData().n_tips > MRCA->getData().n_tips)
                log_supported_by(MRCA, nd);
            else
                for(auto nd2 = MRCA;nd2 and nd2->getData().n_tips == MRCA->getData().n_tips;nd2 = nd2->getParent())
                    log_partial_path_of(nd2, nd);
        }

        conflicts.clear();
        for(auto nd: nodes)
            // If we have (a) some, but not all of the include group
            //            (b) any of the exclude group
            if (n_include_tips(nd) < n_tips(nd) and n_include_tips(nd) < L2)
                conflicts.push_back(nd);

        for(auto nd: nodes)
            n_include_tips(nd) = 0;
            
        for(auto conflicting_node: conflicts)
            log_conflicts_with(conflicting_node, nd);

        if (conflicts.empty() and conflicts_or_resolved_by)
            log_resolved_by(MRCA, nd);

#ifdef CHECK_MARKS
        for(const auto nd2: iter_post_const(*induced_tree2))
            assert(n_include_tips(nd2) == 0);
#endif

        // nd -> MRCA
        if (not conflicts_or_resolved_by)
        {
            summary_node(nd) = MRCA;

            destroy_children(nd);
            destroy_children(MRCA);
        }
    }
}

struct DisplayedStatsState : public TaxonomyDependentTreeProcessor<Tree_t> {
    json document;
    std::unique_ptr<Tree_t> summaryTree;
    map<string,string> monotypic_nodes;
    map<long,const Tree_t::node_type*> taxOttIdToNode;
    map<long,Tree_t::node_type*> summaryOttIdToNode;
    map<long,const Tree_t::node_type*> constSummaryOttIdToNode;
    map<string, Map<string,string>> supported_by;
    map<string, Map<string,string>> partial_path_of;
    map<string, Map<string,string>> conflicts_with;
    map<string, Map<string,string>> resolves;
    map<string, Map<string,string>> resolved_by;
    map<string, Map<string,string>> terminal;
    map<string, set<pair<string, string>>> supported_by_set;
    map<string, set<pair<string, string>>> partial_path_of_set;
    map<string, set<pair<string, string>>> conflicts_with_set;
    map<string, set<pair<string, string>>> resolves_set;
    map<string, set<pair<string, string>>> resolved_by_set;
    map<string, set<pair<string, string>>> terminal_set;
    int numErrors = 0;
    bool treatTaxonomyAsLastTree = false;
    bool headerEmitted = false;

    virtual ~DisplayedStatsState(){}

    void set_terminal(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const string& source)
    {
        add_element(terminal, terminal_set, synth_node, input_node, source);
    }

    void set_supported_by(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const string& source)
    {
        add_element(supported_by, supported_by_set, synth_node, input_node, source);
    }

    void set_partial_path_of(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const string& source)
    {
        add_element(partial_path_of, partial_path_of_set, synth_node, input_node, source);
    }

    void set_conflicts_with(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const string& source)
    {
        add_element(conflicts_with, conflicts_with_set, synth_node, input_node, source);
    }

    void set_resolved_by(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const string& source)
    {
        add_element(resolved_by, resolved_by_set, synth_node, input_node, source);
    }

    void set_resolves(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const string& source)
    {
        add_element(resolves, resolves_set, synth_node, input_node, source);
    }

    bool summarize(OTCLI &otCLI) override {
        if (treatTaxonomyAsLastTree) {
            mapNextTree(otCLI, *taxonomy, true);
        }

        document["num_tips"] = countLeaves(*summaryTree);
        document["root_ott_id"] = summaryTree->getRoot()->getOttId();

        json nodes;
        for(auto nd: iter_post_const(*summaryTree))
        {
            json node;
            string name = nd->getName();
//            if (nd->hasOttId())
//                name = "ott" + std::to_string(nd->getOttId());

            set_support_blob_as_single_element(node, terminal, "terminal", name);
            set_support_blob_as_single_element(node, supported_by, "supported_by", name);
            set_support_blob_as_single_element(node, partial_path_of, "partial_path_of", name);
            set_support_blob_as_array(node, conflicts_with, "conflicts_with", name);
            set_support_blob_as_single_element(node, resolves, "resolves", name);
            set_support_blob_as_array(node, resolved_by, "resolved_by", name);

            if (not node.empty())
                nodes[name] = node;
        }

        // Copy support information to monotypic nodes from their first non-monotypic descendant
        for(const auto& m: monotypic_nodes)
            if (nodes.find(m.second) != nodes.end())
                nodes[m.first] = nodes[m.second];

        document["nodes"] = nodes;
        std::cout<<document.dump(1)<<std::endl;
        return true;
    }

    void mapNextTree(OTCLI & , const Tree_t & tree, bool ) //isTaxoComp is third param
    {
        string source_name = source_from_tree_name(tree.getName());
        document["sources"].push_back(source_name);


        auto ottid_to_node = get_ottid_to_const_node_map(tree);

        {
            auto log_supported_by    = [this, &source_name](const node_t* node2, const node_t* node1) {set_supported_by(node2,node1,source_name);};
            auto log_partial_path_of = [this, &source_name](const node_t* node2, const node_t* node1) {set_partial_path_of(node2,node1,source_name);};
            auto log_conflicts_with  = [this, &source_name](const node_t* node2, const node_t* node1) {set_conflicts_with(node2,node1,source_name);};
            auto log_resolved_by     = [this, &source_name](const node_t* node2, const node_t* node1) {set_resolved_by(node2,node1,source_name);};
            auto log_terminal        = [this, &source_name](const node_t* node2, const node_t* node1) {set_terminal(node2,node1,source_name);};

            perform_conflict_analysis(tree, ottid_to_node,
                                      *summaryTree, constSummaryOttIdToNode,
                                      log_supported_by,
                                      log_partial_path_of,
                                      log_conflicts_with,
                                      log_resolved_by,
                                      log_terminal);
        }
        {
            auto log_supported_by    = [this, &source_name](const node_t* node2, const node_t* node1) {};
            auto log_partial_path_of = [this, &source_name](const node_t* node2, const node_t* node1) {};
            auto log_conflicts_with  = [this, &source_name](const node_t* node2, const node_t* node1) {};
            auto log_resolved_by     = [this, &source_name](const node_t* node2, const node_t* node1) {set_resolves(node1,node2,source_name);};
            auto log_terminal        = [this, &source_name](const node_t* node2, const node_t* node1) {};

            perform_conflict_analysis(*summaryTree, constSummaryOttIdToNode,
                                      tree, ottid_to_node,
                                      log_supported_by,
                                      log_partial_path_of,
                                      log_conflicts_with,
                                      log_resolved_by,
                                      log_terminal);
        }
    }

    virtual bool processTaxonomyTree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<Tree_t>::processTaxonomyTree(otCLI);
        otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
        otCLI.getParsingRules().requireOttIds = false;
        for(auto nd: iter_post_const(*taxonomy))
            if (nd->hasOttId())
                taxOttIdToNode[nd->getOttId()] = nd;
        return true;
    }

    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) override {
        assert(taxonomy != nullptr);
        if (summaryTree == nullptr) {
            summaryTree = std::move(tree);
            monotypic_nodes = suppressAndRecordMonotypic(*summaryTree);
            for(auto nd: iter_post(*summaryTree))
                if (nd->hasOttId())
                    summaryOttIdToNode[nd->getOttId()] = nd;
            constSummaryOttIdToNode = get_ottid_to_const_node_map(*summaryTree);
            computeDepth(*summaryTree);
            return true;
        }
        requireTipsToBeMappedToTerminalTaxa(*tree, taxOttIdToNode);
        computeDepth(*tree);
        computeSummaryLeaves(*tree, summaryOttIdToNode);

        mapNextTree(otCLI, *tree, false);
        return true;
    }

};
bool handleCountTaxonomy(OTCLI & otCLI, const std::string &);

bool handleCountTaxonomy(OTCLI & otCLI, const std::string &) {
    DisplayedStatsState * proc = static_cast<DisplayedStatsState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->treatTaxonomyAsLastTree = true;
    return true;
}

int main(int argc, char *argv[]) {
    std::string explanation{"takes at least 2 newick file paths: a taxonomy, a full supertree, and some number of input trees.\n"};
    OTCLI otCLI("otc-annotate-synth",
                explanation.c_str(),
                "taxonomy.tre synth.tre inp1.tre inp2.tre ...");
    otCLI.addFlag('x',
                  "Automatically treat the taxonomy as an input in terms of supporting groups",
                  handleCountTaxonomy,
                  false);
    DisplayedStatsState proc;
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
