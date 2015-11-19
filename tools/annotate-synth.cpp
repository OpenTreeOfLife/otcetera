#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include <tuple>
#include <sstream>
#include <cstring>
#include <unordered_map>

using namespace otc;

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

string getNodeName(const string& name)
{
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

string quote(const string& s)
{
    return '"'+s+'"';
};


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

void trace_clean_marks(Tree_t::node_type* node)
{
    do
    {
        mark(node) = 0;
        node = node->getParent();
    } while (node);
}

void trace_clean_marks_from_synth(const Tree_t& tree)
{
    for(auto leaf: iter_leaf_const(tree))
    {
        auto leaf2 = summary_node(leaf);
        trace_clean_marks(leaf2);
    }
}

string source_node_name(const Tree_t::node_type* input_node, const Tree_t& input_tree)
{
    string source = quote(source_from_tree_name(input_tree.getName()));
    string node_in_study = quote(getNodeName(input_node->getName()));
    std::ostringstream study_tree_node;;
    study_tree_node<<"["<<source<<", "<<node_in_study<<"]";
    return study_tree_node.str();
}

struct DisplayedStatsState : public TaxonomyDependentTreeProcessor<Tree_t> {
    std::unique_ptr<Tree_t> summaryTree;
    std::map<long,const Tree_t::node_type*> taxOttIdToNode;
    std::map<long,Tree_t::node_type*> summaryOttIdToNode;
    std::unordered_multimap<string,string> supported_by;
    std::unordered_multimap<string,string> partial_path_of;
    std::unordered_multimap<string,string> conflicts_with;
    std::unordered_multimap<string,string> could_resolve;
    std::unordered_multimap<string,string> terminal;
    int numErrors = 0;
    bool treatTaxonomyAsLastTree = false;
    bool headerEmitted = false;
    int numTrees = 0;
    virtual ~DisplayedStatsState(){}

    void set_terminal(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        terminal.insert({synth_node->getName(), source_node_name(input_node,input_tree)});
    }

    void set_supported_by(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        supported_by.insert({synth_node->getName(), source_node_name(input_node,input_tree)});
    }

    void set_partial_path_of(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        partial_path_of.insert({synth_node->getName(), source_node_name(input_node,input_tree)});
    }

    void set_conflicts_with(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        conflicts_with.insert({synth_node->getName(), source_node_name(input_node,input_tree)});
    }

    void set_could_resolve(const Tree_t::node_type* synth_node, const Tree_t::node_type* input_node, const Tree_t& input_tree)
    {
        could_resolve.insert({synth_node->getName(), source_node_name(input_node,input_tree)});
    }

    bool summarize(OTCLI &otCLI) override {
        if (treatTaxonomyAsLastTree) {
            mapNextTree(otCLI, *taxonomy, true);
        }

        std::cout<<"  \"nodes\": {\n";
        for(auto nd: iter_post_const(*summaryTree))
        {
            string name = nd->getName();
            if (nd->hasOttId())
                name = "ott" + std::to_string(nd->getOttId());

            int sc = supported_by.count(name);
            int tc = terminal.count(name);
            int pc = partial_path_of.count(name);
            int cc = conflicts_with.count(name);
            if (sc + tc + pc + cc == 0) continue;

            std::cout<<"    "<<quote(name)<<": { \n";
            if (sc)
            {
                std::cout<<"      \"supported-by\": [ ";
                auto range = supported_by.equal_range(name);
                auto start = range.first;
                auto end   = range.second;
                for (auto iter = start; iter!=end; ++iter)
                {
                    if (iter != start) std::cout<<"                        ";
                    std::cout<<iter->second;
                    auto next = iter; ++next;
                    if (next != end)
                        std::cout<<",\n";
                }
                std::cout<<" ]";
                if (tc + pc + cc > 0)
                    std::cout<<",";
                std::cout<<"\n";
            }
            if (tc)
            {
                std::cout<<"      \"terminal\": [ ";
                auto range = terminal.equal_range(name);
                auto start = range.first;
                auto end   = range.second;
                for (auto iter = start; iter!=end; ++iter)
                {
                    if (iter != start) std::cout<<"                        ";
                    std::cout<<iter->second;
                    auto next = iter; ++next;
                    if (next != end)
                        std::cout<<",\n";
                }
                std::cout<<" ]";
                if (pc + cc > 0)
                    std::cout<<",";
                std::cout<<"\n";
            }
            if (pc)
            {
                std::cout<<"      \"partial-path-of\": [ ";
                auto range = partial_path_of.equal_range(name);
                auto start = range.first;
                auto end   = range.second;
                for (auto iter = start; iter!=end; ++iter)
                {
                    if (iter != start) std::cout<<"                        ";
                    std::cout<<iter->second;
                    auto next = iter; ++next;
                    if (next != end)
                        std::cout<<",\n";
                }
                std::cout<<" ]";
                if (cc > 0)
                    std::cout<<",";
                std::cout<<"\n";
            }
            if (cc)
            {
                std::cout<<"      \"conflicts-with\": [ ";
                auto range = conflicts_with.equal_range(name);
                auto start = range.first;
                auto end   = range.second;
                for (auto iter = start; iter!=end; ++iter)
                {
                    if (iter != start) std::cout<<"                        ";
                    std::cout<<iter->second;
                    auto next = iter; ++next;
                    if (next != end)
                        std::cout<<",\n";
                }
                std::cout<<" ]\n";
            }
            std::cout<<"    }\n";
        }
        std::cout<<"  }\n";

        return true;
    }

    void mapNextTree(OTCLI & otCLI, const Tree_t & tree, bool isTaxoComp)
    {
        for(const auto nd: iter_post_const(tree))
        {
            if (not nd->getParent()) continue;

            auto MRCA_include = trace_include_group_find_MRCA(nd, 1);
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
                if (MRCA_include->getParent() and mark(MRCA_include->getParent()) == 2)
                    set_supported_by(MRCA_include, nd, tree);
                else
                {
                    for(auto path_node = MRCA_include;path_node and not is_marked(path_node,2);path_node = path_node->getParent())
                        set_partial_path_of(path_node, nd, tree);
                }
            }

            trace_clean_marks_from_synth(tree);
#ifndef CHECK_MARKS
            for(const auto nd2: iter_post_const(*summaryTree))
            {
                assert(mark(nd2) == 0);
            }
#endif
        }
        numTrees += 1;
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
