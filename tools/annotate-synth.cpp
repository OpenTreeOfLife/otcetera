// See https://github.com/OpenTreeOfLife/opentree/wiki/Open-Tree-of-Life-APIs-v3#synthetic-tree
#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include <tuple>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include "json.hpp"

using namespace otc;
using json = nlohmann::json;

using std::vector;
using std::string;
using std::map;

struct RTNodeDepth
{
    int depth = 0;
    int mark = 0;
    RootedTreeNode<RTNodeDepth>* summary_node;
};

using Tree_t = RootedTree<RTNodeDepth, RTreeNoData>;
using node_t = Tree_t::node_type;


int depth(const Tree_t::node_type* node)
{
    return node->getData().depth;
}

int& depth(Tree_t::node_type* node)
{
    return node->getData().depth;
}

int mark(const Tree_t::node_type* node)
{
    return node->getData().mark;
}

int& mark(Tree_t::node_type* node)
{
    return node->getData().mark;
}

bool is_marked(const Tree_t::node_type* node, int bits)
{
    return (mark(node)&bits) == bits;
}

void set_mark(Tree_t::node_type* node, int bits)
{
    mark(node) |= bits;
}

Tree_t::node_type* summary_node(const Tree_t::node_type* node)
{
    return node->getData().summary_node;
}

Tree_t::node_type*& summary_node(Tree_t::node_type* node)
{
    return node->getData().summary_node;
}

Tree_t::node_type* nmParent(Tree_t::node_type* node)
{
    do
    {
        node = node->getParent();
    } while (node and node->isOutDegreeOneNode());
    return node;
}

void computeDepth(Tree_t& tree)
{
    tree.getRoot()->getData().depth = 1;
    for(auto nd: iter_pre(tree))
    {
        if (not nd->getParent()) continue;

        nd->getData().depth = nd->getParent()->getData().depth + 1;
    }
}

void computeSummaryNodes(Tree_t& tree, const map<long,Tree_t::node_type*>& summaryOttIdToNode)
{
    for(auto leaf: iter_leaf(tree))
        summary_node(leaf) = summaryOttIdToNode.at(leaf->getOttId());
}

long n_leaves(const node_t* nd)
{
    long count = 0;
    for(auto nd2: iter_post_n_const(*nd))
        if (nd2->isTip())
            count++;
    return count;
}

long n_leaves(const Tree_t& t)
{
    return n_leaves(t.getRoot());
}

long n_children(const node_t* nd)
{
    long count = 0;
    for(auto nd2: iter_child_const(*nd))
        count++;
    return count;
}

long n_marked_children(const node_t* nd, int bits)
{
    long count = 0;
    for(auto nd2: iter_child_const(*nd))
        if (is_marked(nd2, bits))
            count++;
    return count;
}

string source_from_tree_name(const string& name)
{
    const char* start = strrchr(name.c_str(),' ');
    start++;
    assert(start);

    const char* end = strrchr(name.c_str(),'.');
    return name.substr(start - name.c_str(), end-start);
}

string study_from_tree_name(const string& name)
{
    const char* start = strrchr(name.c_str(),' ');
    start++;
    assert(start);

    const char* end = strrchr(name.c_str(),'_');
    return name.substr(start - name.c_str(), end-start);
}


string tree_in_study_from_tree_name(const string& name)
{
    const char* start = strrchr(name.c_str(),'_');
    start++;
    assert(start);

    const char* end = strrchr(name.c_str(),'.');
    return name.substr(start - name.c_str(), end-start);
}

string getNodeName(const Tree_t::node_type* node)
{
    string name = node->getName();
    if (node->hasOttId())
        name = "ott" + std::to_string(node->getOttId());

    const char* start = name.c_str();
    const char* end = name.c_str() + name.size();
    while(strchr(" \t_",*start) and start < end)
        start++;
    while(strchr(" \t_",*(end-1)) and start < end)
        end--;
    if (start == end)
        throw OTCError()<<"Node name '"<<name<<"' contracted to nothing!";
    return name.substr(start - name.c_str(), end-start);
}

Tree_t::node_type* trace_to_parent(Tree_t::node_type* node, int bits)
{
    // node should be already marked
    assert(is_marked(node, bits));
    
    // move to parent and mark it.
    node = node->getParent();
    set_mark(node, bits);

    return node;
}

Tree_t::node_type* get_root(Tree_t::node_type* node)
{
    while(node->getParent())
        node = node->getParent();
    return node;
}

const Tree_t::node_type* get_root(const Tree_t::node_type* node)
{
    while(node->getParent())
        node = node->getParent();
    return node;
}

/// Walk up the tree from node1 and node2 until we find the common ancestor, marking all the way.
Tree_t::node_type* trace_find_MRCA(Tree_t::node_type* node1, Tree_t::node_type* node2, int bits1, int bits2)
{
    assert(node1 or node2);
    if (not node1)
    {
        assert(is_marked(node2, bits2));
        return node2;
    }

    if (not node2)
    {
        assert(is_marked(node1, bits1));
        return node1;
    }

    assert(node1 and node2);
    assert(get_root(node1) == get_root(node2));

    assert(is_marked(node1, bits1));
    assert(is_marked(node2, bits2));
    
    while(depth(node1) > depth(node2))
        node1 = trace_to_parent(node1, bits1);
    while(depth(node1) < depth(node2))
        node2 = trace_to_parent(node2, bits2);

    assert(depth(node1) == depth(node2));
    while(node1 != node2)
    {
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

Tree_t::node_type* trace_include_group_find_MRCA(const Tree_t::node_type* node, int bits)
{
    Tree_t::node_type* MRCA = nullptr;
    for(auto leaf: iter_leaf_n_const(*node))
    {
        auto leaf2 = summary_node(leaf);
        mark(leaf2) |= bits;
        MRCA = trace_find_MRCA(MRCA, leaf2, bits, bits);
        assert(is_marked(MRCA,bits));
        if (MRCA->hasChildren())
            assert(n_marked_children(MRCA,bits)>0);
    }
    return MRCA;
}


Tree_t::node_type* trace_exclude_group_find_MRCA(const Tree_t::node_type* node, int bits1, int bits2)
{
    // Using bits1 to exclude leafs seems like a dumb way to iterate over excluded leaves.
    node = get_root(node);

    Tree_t::node_type* MRCA = nullptr;
    for(auto leaf: iter_leaf_n(*node))
    {
        auto leaf2 = summary_node(leaf);

        if (is_marked(leaf2, bits1)) continue;

        mark(leaf2) |= bits2;
        MRCA = trace_find_MRCA(MRCA, leaf2, bits2, bits2);
    }
    return MRCA;
}

void find_conflicts(Tree_t::node_type* node, vector<Tree_t::node_type*>& conflicts)
{
    while (node and mark(node))
    {
        if (is_marked(node,1) and is_marked(node,2) and node->getParent() and is_marked(node->getParent(),1))
            conflicts.push_back(node);
        node = node->getParent();
    } 
}

void find_conflicts(const Tree_t& tree, vector<Tree_t::node_type*>& conflicts)
{
    conflicts.clear();
    for(auto leaf: iter_leaf_const(tree))
    {
        auto leaf2 = summary_node(leaf);
        find_conflicts(leaf2, conflicts);
    }
}

void trace_clean_marks(Tree_t::node_type* node)
{
    while (node and mark(node))
    {
        mark(node) = 0;
        node = node->getParent();
    } 
}

void trace_clean_marks_from_synth(const Tree_t& tree)
{
    for(auto leaf: iter_leaf_const(tree))
    {
        auto leaf2 = summary_node(leaf);
        trace_clean_marks(leaf2);
    }
}

json source_node(const Tree_t::node_type* input_node, const Tree_t& input_tree)
{
    string source = source_from_tree_name(input_tree.getName());
    string node_in_study = getNodeName(input_node);
    return {source,node_in_study};
}

struct DisplayedStatsState : public TaxonomyDependentTreeProcessor<Tree_t> {
    json document;
    std::unique_ptr<Tree_t> summaryTree;
    std::map<long,const Tree_t::node_type*> taxOttIdToNode;
    std::map<long,Tree_t::node_type*> summaryOttIdToNode;
    std::unordered_multimap<string,json> supported_by;
    std::unordered_multimap<string,json> partial_path_of;
    std::unordered_multimap<string,json> conflicts_with;
    std::unordered_multimap<string,json> could_resolve;
    std::unordered_multimap<string,json> terminal;
    int numErrors = 0;
    bool treatTaxonomyAsLastTree = false;
    bool headerEmitted = false;

    virtual ~DisplayedStatsState(){}

    void set_terminal(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        terminal.insert({synth_node->getName(), source_node(input_node,input_tree)});
    }

    void set_supported_by(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        supported_by.insert({synth_node->getName(), source_node(input_node,input_tree)});
    }

    void set_partial_path_of(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        partial_path_of.insert({synth_node->getName(), source_node(input_node,input_tree)});
    }

    void set_conflicts_with(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        conflicts_with.insert({synth_node->getName(), source_node(input_node,input_tree)});
    }

    void set_could_resolve(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        could_resolve.insert({synth_node->getName(), source_node(input_node,input_tree)});
    }

    bool summarize(OTCLI &otCLI) override {
        if (treatTaxonomyAsLastTree) {
            mapNextTree(otCLI, *taxonomy, true);
        }

//        document["date_completed"] = "date";
//        document["tree_id"] = "an idenftifier";
//        document["taxonomy_version"] = "2.9draft12";
        document["num_tips"] = n_leaves(*summaryTree);
//        document["run_time"] = "an estimate of the time taken to build the tree";
//        document["num_source_trees"] = numTrees;
//        document["num_source_studies"] = studies.size();
//        document["root_taxon_name"] = "life";
        document["root_ott_id"] = summaryTree->getRoot()->getOttId();
//        document["generated_by"] = "propinquity";
//        document["filtered_flags"] = "list of taxon flags";

        json nodes;
        for(auto nd: iter_post_const(*summaryTree))
        {
            json node;
            string name = nd->getName();
            if (nd->hasOttId())
                name = "ott" + std::to_string(nd->getOttId());

            {
                json j_supported_by = json::array();

                auto range = supported_by.equal_range(name);
                for (auto iter = range.first; iter!=range.second; ++iter)
                    j_supported_by.push_back(iter->second);
                if (not j_supported_by.empty())
                    node["supported_by"] = j_supported_by;
            }
            {
                json j_terminal = json::array();

                auto range = terminal.equal_range(name);
                for (auto iter = range.first; iter!=range.second; ++iter)
                    j_terminal.push_back(iter->second);
                if (not j_terminal.empty())
                    node["terminal"] = j_terminal;
            }
            {
                json j_partial_path_of = json::array();

                auto range = partial_path_of.equal_range(name);
                for (auto iter = range.first; iter!=range.second; ++iter)
                    j_partial_path_of.push_back(iter->second);
                if (not j_partial_path_of.empty())
                    node["partial_path_of"] = j_partial_path_of;
            }
            {
                json j_conflicts_with = json::array();

                auto range = conflicts_with.equal_range(name);
                for (auto iter = range.first; iter!=range.second; ++iter)
                    j_conflicts_with.push_back(iter->second);
                if (not j_conflicts_with.empty())
                    node["conflicts_with"] = j_conflicts_with;
            }
            {
                json j_could_resolve = json::array();

                auto range = could_resolve.equal_range(name);
                for (auto iter = range.first; iter!=range.second; ++iter)
                    j_could_resolve.push_back(iter->second);
                if (not j_could_resolve.empty())
                    node["could_resolve"] = j_could_resolve;
            }
            if (not node.empty())
                nodes[name] = node;
        }
        document["nodes"] = nodes;
        std::cout<<document.dump(4)<<std::endl;
        return true;
    }

    void mapNextTree(OTCLI & otCLI, const Tree_t & tree, bool isTaxoComp)
    {
        vector<Tree_t::node_type*> conflicts;
        string source_name = source_from_tree_name(tree.getName());
        document["sources"].push_back(source_name);
        
        for(const auto nd: iter_post_const(tree))
        {
            if (not nd->getParent()) continue;

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
        requireTipsToBeMappedToTerminalTaxa(*tree, taxOttIdToNode);
        computeDepth(*tree);
        computeSummaryNodes(*tree, summaryOttIdToNode);

        mapNextTree(otCLI, *tree, false);
        return true;
    }

};


bool handleCountTaxonomy(OTCLI & otCLI, const std::string &) {
    DisplayedStatsState * proc = static_cast<DisplayedStatsState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->treatTaxonomyAsLastTree = true;
    return true;
}

int main(int argc, char *argv[]) {
    std::string explanation{"takes at least 2 newick file paths: a taxonomy,  a full supertree, and some number of input trees.\n"};
    OTCLI otCLI("otc-displayed-stats",
                explanation.c_str(),
                "taxonomy.tre synth.tre inp1.tre inp2.tre ...");
    DisplayedStatsState proc;
    otCLI.addFlag('x',
                  "Automatically treat the taxonomy as an input in terms of supporting groups",
                  handleCountTaxonomy,
                  false);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
