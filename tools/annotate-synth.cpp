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
using Set = std::unordered_set<T>;
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
    int mark = 0;
    RootedTreeNode<RTNodeDepth>* summary_node;
};

using Tree_t = RootedTree<RTNodeDepth, RTreeNoData>;
using node_t = Tree_t::node_type;

int depth(const Tree_t::node_type* node);
int& depth(Tree_t::node_type* node);
Tree_t::node_type* summary_node(const Tree_t::node_type* node);
Tree_t::node_type*& summary_node(Tree_t::node_type* node);
Tree_t::node_type* nmParent(Tree_t::node_type* node);
void computeSummaryLeaves(Tree_t& tree, const map<long,Tree_t::node_type*>& summaryOttIdToNode);
string getSourceNodeNameIfAvailable(const Tree_t::node_type* node);
Tree_t::node_type* trace_to_parent(Tree_t::node_type* node, int bits);
Tree_t::node_type* get_root(Tree_t::node_type* node);
const Tree_t::node_type* get_root(const Tree_t::node_type* node);
Tree_t::node_type* trace_find_MRCA(Tree_t::node_type* node1, Tree_t::node_type* node2, int bits1, int bits2);
Tree_t::node_type* trace_include_group_find_MRCA(const Tree_t::node_type* node, int bits);
Tree_t::node_type* trace_exclude_group_find_MRCA(const Tree_t::node_type* node, int bits1, int bits2);
void find_anc_conflicts(Tree_t::node_type* node, vector<Tree_t::node_type*>& conflicts);
void find_conflicts(const Tree_t& tree, vector<Tree_t::node_type*>& conflicts);
void trace_clean_marks(Tree_t::node_type* node);
void trace_clean_marks_from_synth(const Tree_t& tree);

bool prune_unrecognized = true;
bool handlePruneUnrecognizedTips(OTCLI & otCLI, const std::string & arg) {
    prune_unrecognized = false;
    return true;
}

inline int depth(const Tree_t::node_type* node) {
    return node->getData().depth;
}


inline int& depth(Tree_t::node_type* node) {
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

// uses the OTT Ids in `tree` to fill in the `summary_node` field of each leaf
void computeSummaryLeaves(Tree_t& tree, const map<long,Tree_t::node_type*>& summaryOttIdToNode) {
    for(auto leaf: iter_leaf(tree)) {
        summary_node(leaf) = summaryOttIdToNode.at(leaf->getOttId());
    }
}


string getSourceNodeNameIfAvailable(const Tree_t::node_type* node) {
    string name = node->getName();
    try
    {
        return getSourceNodeName(name);
    }
    catch (std::exception& e)
    {
        return name;
    }
}

// assumes `node` is marked with `bits`. Returns the parent of `node` and
//  also marks the parent with `bits` as a side effect.
Tree_t::node_type* trace_to_parent(Tree_t::node_type* node, int bits) {
    assert(is_marked(node, bits));
    node = node->getParent(); // move to parent and mark it.
    set_mark(node, bits);
    return node;
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

/// Walk up the tree from node1 and node2 until we find the common ancestor, marking all the way.
// for X=1,2 and Y = 3-X
//      assumes nodeX is marked by bitsX 
//      returns the MRCA, and guarantees that the path from MRCA to nodeX is marked with bitsX
//  if nodeX is nullptr, returns the other nodeY marked by bitsY
Tree_t::node_type* trace_find_MRCA(Tree_t::node_type* node1, Tree_t::node_type* node2, int bits1, int bits2) {
    assert(node1 or node2);
    if (not node1) {
        assert(is_marked(node2, bits2));
        return node2;
    }
    if (not node2) {
        assert(is_marked(node1, bits1));
        return node1;
    }
    assert(node1 and node2);
    assert(get_root(node1) == get_root(node2));
    assert(is_marked(node1, bits1));
    assert(is_marked(node2, bits2));
    while (depth(node1) > depth(node2)){
        node1 = trace_to_parent(node1, bits1);
    }
    while (depth(node1) < depth(node2)){
        node2 = trace_to_parent(node2, bits2);
    }
    assert(depth(node1) == depth(node2));
    while (node1 != node2) {
        assert(node1->getParent());
        assert(node2->getParent());
        node1 = trace_to_parent(node1, bits1);
        node2 = trace_to_parent(node2, bits2);
    }
    assert(node1 == node2);
    assert(is_marked(node1, bits1));
    assert(is_marked(node2, bits2));
    return node1;
}

// traces the summary nodes that correspond to the leaves of `node` back to their MRCA
//  on the summary tree. All of the nodes in this induced tree will be flagged by setting
//  `bits` to 1.
// returns the MRCA node (in the summary tree)
Tree_t::node_type* trace_include_group_find_MRCA(const Tree_t::node_type* node, int bits) {
    Tree_t::node_type* MRCA = nullptr;
    for (auto leaf: iter_leaf_n_const(*node)) {
        auto leaf2 = summary_node(leaf);
        mark(leaf2) |= bits;
        MRCA = trace_find_MRCA(MRCA, leaf2, bits, bits);
        assert(is_marked(MRCA, bits));
        if (MRCA->hasChildren()) {
            assert(countMarkedChildren(MRCA,bits)>0);
        }
    }
    return MRCA;
}


// Assumes that the version of the summary tree induced by the include group of `node` has
//  already been marked with `bits1`.
// Returns the MRCA of the exclude group, and assures that every node in the version of the
//  summary tree induced by the exclude group is flagged with `bits2` 
Tree_t::node_type* trace_exclude_group_find_MRCA(const Tree_t::node_type* node, int bits1, int bits2) {
    // Using bits1 to exclude leafs seems like a dumb way to iterate over excluded leaves.
    node = get_root(node);
    Tree_t::node_type* MRCA = nullptr;
    for(auto leaf: iter_leaf_n(*node)) {
        auto leaf2 = summary_node(leaf);
        if (!is_marked(leaf2, bits1)) {
            mark(leaf2) |= bits2;
            MRCA = trace_find_MRCA(MRCA, leaf2, bits2, bits2);
        }
    }
    return MRCA;
}

// Walks from `node` to the induced root (real root or the deepest node with some mark).
// adds the a node to `conflicts` if that node is marked by the include AND exclude bit
//  but the node is NOT the deepest node marked with the include flag.
// The latter (but...NOT) condition avoids flagging nodes that are polytomies which
//  could be resolved in favor of a node as "conflict" nodes.
// This is sub-optimal, because we walk each branch n times if it has n children.
// There a conflict will be discovered n times if it has n children.
// We currently hack around this by checking for duplicate entries in the hash.
void find_anc_conflicts(Tree_t::node_type* node, vector<Tree_t::node_type*>& conflicts) {
    while (node and mark(node)) {
        if (is_marked(node,1)
            and is_marked(node,2)
            and node->getParent()
            and is_marked(node->getParent(), 1)) {
            conflicts.push_back(node);
        }
        node = node->getParent();
    } 
}

// calls find_anc_conflicts for each leaf to fill `conflicts`
void find_conflicts(const Tree_t& tree, vector<Tree_t::node_type*>& conflicts) {
    conflicts.clear();
    for(auto leaf: iter_leaf_const(tree)) {
        auto leaf2 = summary_node(leaf);
        find_anc_conflicts(leaf2, conflicts);
    }
}


void trace_clean_marks(Tree_t::node_type* node) {
    while (node and mark(node)) {
        mark(node) = 0;
        node = node->getParent();
        // MTH should  be able to bail out here if the next node is already mark(node) == 0, right?
    } 
}

// sets marks to 0 for all nodes of the summary tree that are induced by the leaf set of tree.
void trace_clean_marks_from_synth(const Tree_t& tree) {
    for(auto leaf: iter_leaf_const(tree)) {
        auto leaf2 = summary_node(leaf);
        trace_clean_marks(leaf2);
    }
}

void pruneUnmapped(Tree_t& tree, const map<long,const Tree_t::node_type*>& taxOttIdToNode)
{
    vector<Tree_t::node_type*> remove;
    for(auto nd: iter_leaf(tree))
    {
        long id = nd->getOttId();
        if (not taxOttIdToNode.count(id))
            remove.push_back(nd);
    }

    for(auto nd: remove)
    {
        auto parent = nd->getParent();
        pruneAndDelete(tree, nd);
        while(parent and not parent->hasOttId() and parent->isTip())
        {
            auto tmp = parent;
            parent = parent->getParent();
            pruneAndDelete(tree, tmp);
        }
    }
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

void add_element(map<string, Map<string, string>>& m, map<string, Set<pair<string,string>>>& s,
                 const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
{
    string synth = synth_node->getName();
    string source = source_from_tree_name(input_tree.getName());
    string node = getSourceNodeNameIfAvailable(input_node);

    pair<string,string> x{source, node};

    if (not s[synth].count(x))
    {
        s[synth].insert(x);
        m[synth].insert(x);
    }
}

struct DisplayedStatsState : public TaxonomyDependentTreeProcessor<Tree_t> {
    json document;
    std::unique_ptr<Tree_t> summaryTree;
    map<long,const Tree_t::node_type*> taxOttIdToNode;
    map<long,Tree_t::node_type*> summaryOttIdToNode;
    map<string, Map<string,string>> supported_by;
    map<string, Map<string,string>> partial_path_of;
    map<string, Map<string,string>> conflicts_with;
    map<string, Map<string,string>> could_resolve;
    map<string, Map<string,string>> terminal;
    map<string, Set<pair<string, string>>> supported_by_set;
    map<string, Set<pair<string, string>>> partial_path_of_set;
    map<string, Set<pair<string, string>>> conflicts_with_set;
    map<string, Set<pair<string, string>>> could_resolve_set;
    map<string, Set<pair<string, string>>> terminal_set;
    int numErrors = 0;
    bool treatTaxonomyAsLastTree = false;
    bool headerEmitted = false;

    virtual ~DisplayedStatsState(){}

    void set_terminal(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        add_element(terminal, terminal_set, synth_node, input_node, input_tree);
    }

    void set_supported_by(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        add_element(supported_by, supported_by_set, synth_node, input_node, input_tree);
    }

    void set_partial_path_of(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        add_element(partial_path_of, partial_path_of_set, synth_node, input_node, input_tree);
    }

    void set_conflicts_with(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        add_element(conflicts_with, conflicts_with_set, synth_node, input_node, input_tree);
    }

    void set_could_resolve(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        add_element(could_resolve, could_resolve_set, synth_node, input_node, input_tree);
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
            set_support_blob_as_array(node, partial_path_of, "partial_path_of", name);
            set_support_blob_as_array(node, conflicts_with, "conflicts_with", name);
            set_support_blob_as_array(node, could_resolve, "could_resolve", name);

            if (not node.empty())
                nodes[name] = node;
        }
        document["nodes"] = nodes;
        std::cout<<document.dump(4)<<std::endl;
        return true;
    }

    void mapNextTree(OTCLI & , const Tree_t & tree, bool ) //isTaxoComp is third param
    {
        vector<Tree_t::node_type*> conflicts;
        string source_name = source_from_tree_name(tree.getName());
        document["sources"].push_back(source_name);
        
        for(const auto nd: iter_post_const(tree))
        {
            if (not nd->getParent()) continue;

            // Ignore knuckles in input trees.
            //
            // Note that in general, if we've pruned this tree down to match the shared taxon set
            // then this could produce knuckles.
            if (nd->isOutDegreeOneNode()) continue;

            auto MRCA_include = trace_include_group_find_MRCA(nd, 1);
            assert(not is_marked(MRCA_include,2));
            auto MRCA_exclude = trace_exclude_group_find_MRCA(nd, 1, 2);

            MRCA_exclude = trace_find_MRCA(MRCA_include, MRCA_exclude, 0, 2);

            if (nd->isTip())
            {
                assert(mark(MRCA_include) == 1);
                for(auto path_node = MRCA_include;path_node and not is_marked(path_node,2);path_node = path_node->getParent())
                    set_terminal(path_node, nd, tree);
            }
            else if (mark(MRCA_include) == 1)
            {
                if (MRCA_include->getParent() and mark(nmParent(MRCA_include)) == 2)
                {
                    auto path_node = MRCA_include;
                    do {
                        set_supported_by(path_node, nd, tree);
                        path_node = path_node->getParent();
                    } while(path_node->isOutDegreeOneNode());
                }
                else
                {
                    for(auto path_node = MRCA_include;path_node and not is_marked(path_node,2);path_node = path_node->getParent())
                        set_partial_path_of(path_node, nd, tree);
                }
            }
            assert(is_marked(MRCA_include,1));
            bool conflicts_or_could_resolve = is_marked(MRCA_include,2);

            find_conflicts(tree, conflicts);
            trace_clean_marks_from_synth(tree);
            
            if (nd->isTip() or mark(MRCA_include) == 1) assert(conflicts.empty());

            for(auto conflicting_node: conflicts)
                set_conflicts_with(conflicting_node, nd, tree);

            if (conflicts.empty() and conflicts_or_could_resolve)
                set_could_resolve(MRCA_include, nd, tree);

#ifdef CHECK_MARKS
            for(const auto nd2: iter_post_const(*summaryTree))
            {
                assert(mark(nd2) == 0);
            }
#endif
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
        computeDepth(*tree);
        assert(taxonomy != nullptr);
        if (summaryTree == nullptr) {
            summaryTree = std::move(tree);
            for(auto nd: iter_post(*summaryTree))
                if (nd->hasOttId())
                    summaryOttIdToNode[nd->getOttId()] = nd;
            return true;
        }
        if (prune_unrecognized)
            pruneUnmapped(*tree, taxOttIdToNode);
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
    otCLI.addFlag('p',
                  "Prune input tips that are not part of the supertree.  Defaults to false",
                  handlePruneUnrecognizedTips,
                  true);
    otCLI.addFlag('x',
                  "Automatically treat the taxonomy as an input in terms of supporting groups",
                  handleCountTaxonomy,
                  false);
    DisplayedStatsState proc;
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
