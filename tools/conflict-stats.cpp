// See https://github.com/OpenTreeOfLife/opentree/wiki/Open-Tree-of-Life-APIs-v3#synthetic-tree
#include <iostream>
#include <tuple>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include "otc/conflict.h"
#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"

using namespace otc;
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
using std::unique_ptr;
using std::string;

namespace po = boost::program_options;
using po::variables_map;

// TODO: Could we exemplify tips here, if we had access to the taxonomy?
variables_map parse_cmd_line(int argc,char* argv[]) 
{ 
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("synth", value<string>(),"Filename for the synthesis tree")
        ("input", value<vector<string>>()->composing(),"Filename for input trees")
        ;

    options_description reporting("Reporting options");
    reporting.add_options()
	("each",value<bool>()->default_value(false),"Show separate results for each input tree")
	("all",value<bool>()->default_value(true),"Show accumulated over all input trees")
	("switch","Count synth nodes instead of input tree nodes")
        ("names,N","Write out node names instead of counts.")
	;

    options_description visible;
    visible.add(reporting).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("synth", 1);
    p.add("input", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-conflict-stats <synth-tree> <input tree1> <input tree2> ... [OPTIONS]\n"
                                                    "Count or report input edges of each annotation type.\n\n"
						    "For each relation of the form (synth_node RELATION input_node), find the input nodes\n"
						    " that satisfy it for at least 1 synth_node.\n"
						    "If --switch is given, find the synth_nodes that satisfy it for at least 1 input_node.",
                                                    visible, invisible, p);

    return vm;
}

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

Tree_t::node_type* summary_node(const Tree_t::node_type* node);
Tree_t::node_type*& summary_node(Tree_t::node_type* node);
void computeSummaryLeaves(Tree_t& tree, const map<long,Tree_t::node_type*>& summaryOttIdToNode);
string getSourceNodeNameIfAvailable(const Tree_t::node_type* node);

inline Tree_t::node_type* summary_node(const Tree_t::node_type* node) {
    return node->getData().summary_node;
}

inline Tree_t::node_type*& summary_node(Tree_t::node_type* node) {
    return node->getData().summary_node;
}

// uses the OTT Ids in `tree` to fill in the `summary_node` field of each leaf
void computeSummaryLeaves(Tree_t& tree, const map<long,Tree_t::node_type*>& summaryOttIdToNode) {
    for(auto leaf: iter_leaf(tree)) {
        summary_node(leaf) = summaryOttIdToNode.at(leaf->getOttId());
    }
}

string getSourceNodeNameIfAvailable(const Tree_t::node_type* node) {
    string name = node->getName();
    if (name.empty())
	throw OTCError()<<"Cannot get name for unnamed node!";

    auto source = getSourceNodeName(name);
    if (source)
        return *source;
    else
        return name;
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
        if (nd->getName().size())
        {
            auto child = nd->getFirstChild();
            assert(not child->isOutDegreeOneNode());
            assert(child->getName().size());
            assert(to_child.count(nd->getName()) == 0);
            to_child[nd->getName()] = child->getName();
        }
        delMonotypicNode(nd,tree);
    }
    return to_child;
}

struct stats
{
    set<pair<string,string>> supported_by;
    set<pair<string,string>> partial_path_of;
    set<pair<string,string>> terminal;
    set<pair<string,string>> conflicts_with;
    set<pair<string,string>> resolves;
    set<pair<string,string>> resolved_by;
};

int numErrors = 0;
bool headerEmitted = false;

void add_element(set<pair<string, string>>& s,
                 const Tree_t::node_type*,
		 const Tree_t::node_type* input_node,
		 const string& source)
{
    // We only care about non-monotypic nodes.
    if (input_node->isOutDegreeOneNode()) return;
    
    string node = getSourceNodeNameIfAvailable(input_node);

    s.insert({source, node});
}

void mapNextTree1(const Tree_t& summaryTree,
		  const map<long, const Tree_t::node_type*>& constSummaryOttIdToNode,
		  const Tree_t & tree,
		  stats& s) //isTaxoComp is third param
{
    string source_name = source_from_tree_name(tree.getName());

    auto ottid_to_node = get_ottid_to_const_node_map(tree);
    {
        auto log_supported_by    = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.supported_by,node2,node1,source_name);};
        auto log_partial_path_of = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.partial_path_of,node2,node1,source_name);};
        auto log_conflicts_with  = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.conflicts_with,node2,node1,source_name);};
        auto log_resolved_by     = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.resolved_by,node2,node1,source_name);};
        auto log_terminal        = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.terminal,node2,node1,source_name);};

        perform_conflict_analysis(tree, ottid_to_node,
                                  summaryTree, constSummaryOttIdToNode,
                                  log_supported_by,
                                  log_partial_path_of,
                                  log_conflicts_with,
                                  log_resolved_by,
                                  log_terminal);
    }
    {
        auto nothing    = [](const node_t*, const node_t*) {};
        auto log_resolved_by     = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.resolves,node1,node2,source_name);};

        perform_conflict_analysis(summaryTree, constSummaryOttIdToNode,
                                  tree, ottid_to_node,
                                  nothing,
                                  nothing,
                                  nothing,
                                  log_resolved_by,
                                  nothing);
    }
}

void mapNextTree2(const Tree_t& summaryTree,
		  const map<long, const Tree_t::node_type*>& constSummaryOttIdToNode,
		  const Tree_t & tree,
		  stats& s) //isTaxoComp is third param
{
    string source_name = source_from_tree_name(summaryTree.getName());

    auto ottid_to_node = get_ottid_to_const_node_map(tree);
    {
        auto log_supported_by    = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.supported_by,node2,node1,source_name);};
        auto log_partial_path_of = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.partial_path_of,node2,node1,source_name);};
        auto log_conflicts_with  = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.conflicts_with,node2,node1,source_name);};
        auto log_resolved_by     = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.resolved_by,node2,node1,source_name);};
        auto log_terminal        = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.terminal,node2,node1,source_name);};

        perform_conflict_analysis(tree, ottid_to_node,
                                  summaryTree, constSummaryOttIdToNode,
                                  log_supported_by,
                                  log_partial_path_of,
                                  log_conflicts_with,
                                  log_resolved_by,
                                  log_terminal);
    }
    {
        auto nothing    = [](const node_t*, const node_t*) {};
        auto log_resolved_by     = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.resolves,node1,node2,source_name);};

        perform_conflict_analysis(summaryTree, constSummaryOttIdToNode,
                                  tree, ottid_to_node,
                                  nothing,
                                  nothing,
                                  nothing,
                                  log_resolved_by,
                                  nothing);
    }
}

void mapNextTree(const Tree_t& summaryTree,
		 const map<long, const Tree_t::node_type*>& constSummaryOttIdToNode,
		 const Tree_t & tree,
		 stats& s, //isTaxoComp is third param
		 bool sw)
{
    if (not sw)
	mapNextTree1(summaryTree, constSummaryOttIdToNode, tree, s);
    else
	mapNextTree2(summaryTree, constSummaryOttIdToNode, tree, s);
}

void show_header(std::ostream& o)
{
    o<<"supported_by"<<"\t"<<"partial_path_of"<<"\t"<<"terminal"<<"\t"<<"conflicts_with"<<"\t"<<"resolves"<<"\t"<<"resolved_by"<<"\tname\n";
}

void show_stats(std::ostream& o, const stats& s, const string& name)
{
    o<<s.supported_by.size()<<"\t"<<s.partial_path_of.size()<<"\t"<<s.terminal.size()<<"\t"<<s.conflicts_with.size()<<"\t"<<s.resolves.size()<<"\t"<<s.resolved_by.size()<<"\t"<<name<<"\n";
}

void show_group(std::ostream& o, const set<pair<string,string>>& group, const string& relation, const string& name)
{
    for(auto& x: group)
	o<<"relation = "<<relation<<"   tree = "<<x.first<<"   node = "<<x.second<<"   group = "<<name<<"\n";
}

void show_names(std::ostream& o, const stats& s, const string& name)
{
    show_group(o, s.supported_by, "supported_by", name);
    show_group(o, s.partial_path_of, "partial_path_of", name);
    show_group(o, s.terminal, "terminal", name);
    show_group(o, s.conflicts_with, "conflicts_with", name);
    show_group(o, s.resolves, "resolves", name);
    show_group(o, s.resolved_by, "resolved_by", name);
    return;
}



int main(int argc, char *argv[]) {
    try {
        variables_map args = parse_cmd_line(argc,argv);
        string synthfilename = args["synth"].as<string>();
        vector<string> inputs = args["input"].as<vector<string>>();
	bool each = args["each"].as<bool>();
	bool all = args["all"].as<bool>();
	bool names = args.count("names");
	bool sw = args.count("switch");
	if (not each and not all)
	    throw OTCError()<<"No output requested.";
        
	// 1. Load and process summary tree.
        auto summaryTree = get_tree<Tree_t>(synthfilename);
	auto summaryOttIdToNode = get_ottid_to_node_map(*summaryTree);
	auto constSummaryOttIdToNode = get_ottid_to_const_node_map(*summaryTree);
	auto monotypic_nodes = suppressAndRecordMonotypic(*summaryTree);
	computeDepth(*summaryTree);
	stats global;

	// 2. Load and process input trees.
	if (not names)
	    show_header(std::cout);
        for(const auto& filename: inputs) {
	    auto tree = get_tree<Tree_t>(filename);
	    computeDepth(*tree);

	    computeSummaryLeaves(*tree, summaryOttIdToNode);

	    string source_name = source_from_tree_name(tree->getName());

	    if (all)
		mapNextTree(*summaryTree, constSummaryOttIdToNode, *tree, global, sw);
	    if (each)
	    {
		stats local;
		mapNextTree(*summaryTree, constSummaryOttIdToNode, *tree, local, sw);
		if (names)
		    show_names(std::cout, local, source_name);
		else
		    show_stats(std::cout, local, source_name);
	    }
        }

	// 3. Write out the stats
	if (all and names)
	    show_names(std::cout, global, "ALL");
	else if (all)
	    show_stats(std::cout, global, "ALL");
	
    }
    catch (std::exception& e)
    {
        std::cerr<<"otc-conflict-stats: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}
