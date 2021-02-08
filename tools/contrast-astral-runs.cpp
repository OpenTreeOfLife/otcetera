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
        ("first", value<string>(),"First assembled astral tree file (output followed by inputs in newick).")
        ("second", value<string>(),"Second assembled astral tree file (output followed by inputs in newick).")
        ;

    options_description reporting("Reporting options");
    // reporting.add_options()
    //     ("each",value<bool>()->default_value(false),"Show separate results for each input tree")
    //     ("all",value<bool>()->default_value(true),"Show accumulated over all input trees")
    //     ("switch","Count synth nodes instead of input tree nodes")
    //     ("names,N","Write out node names instead of counts.")
    //     ;

    options_description visible;
    visible.add(reporting).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("first", 1);
    p.add("second", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-contrast-astral-runs <first> <second> [OPTIONS]\n"
                                                    "Calculates conflict stats for astral runs.\n\n"
                                                    "Input files should be newick. Astral output first, followed by the inputs.\n"
                                                    "all trees rooted, and the order of gene trees identical between inputs.",
                                                    visible, invisible, p);

    return vm;
}

namespace std {

template<typename T, typename U>
struct hash<pair<T,U>> {
    std::size_t operator()(const std::pair<T,U>& p) const noexcept {return std::hash<T>()(p.first) * std::hash<U>()(p.second);}
};

} // namespace std

using Tree_t = ConflictTree;
using node_t = Tree_t::node_type;

void compute_summary_leaves(Tree_t& tree, const map<OttId,Tree_t::node_type*>& summaryOttIdToNode);
string get_source_node_name_if_available(const Tree_t::node_type* node);

// uses the OTT Ids in `tree` to fill in the `summary_node` field of each leaf
void compute_summary_leaves(Tree_t& tree, const map<OttId,Tree_t::node_type*>& summaryOttIdToNode) {
    for(auto leaf: iter_leaf(tree)) {
        summary_node(leaf) = summaryOttIdToNode.at(leaf->get_ott_id());
    }
}

string get_source_node_name_if_available(const Tree_t::node_type* node) {
    string name = node->get_name();
    if (name.empty()) {
        throw OTCError() << "Cannot get name for unnamed node!";
    }
    auto source = get_source_node_name(name);
    if (source){
        return *source;
    } else {
        return name;
    }
}

map<string,string> suppress_and_record_monotypic(Tree_t& tree) {
    map<string,string> to_child;
    std::vector<Tree_t::node_type*> remove;
    for (auto nd:iter_post(tree)) {
        if (nd->is_outdegree_one_node()) {
            remove.push_back(nd);
        }
    }

    for (auto nd: remove) {
        if (nd->get_name().size()) {
            auto child = nd->get_first_child();
            assert(not child->is_outdegree_one_node());
            assert(child->get_name().size());
            assert(to_child.count(nd->get_name()) == 0);
            to_child[nd->get_name()] = child->get_name();
        }
        del_monotypic_node(nd,tree);
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
                 const string& source) {
    // We only care about non-monotypic nodes.
    if (input_node->is_outdegree_one_node()) {
        return;
    }
    string node = get_source_node_name_if_available(input_node);
    s.insert({source, node});
}

void mapNextTree1(const Tree_t& summaryTree,
                  const map<OttId, const Tree_t::node_type*>& constSummaryOttIdToNode,
                  const Tree_t & tree,
                  stats& s) {
    typedef Tree_t::node_type node_type;
    std::function<const node_type*(const node_type*,const node_type*)> mrca_of_pair = [](const node_type* n1, const node_type* n2) {return mrca_from_depth(n1,n2);};
    string source_name = source_from_tree_name(tree.get_name());
    auto ottid_to_node = get_ottid_to_const_node_map(tree);
    {
        auto log_supported_by    = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.supported_by,node2,node1,source_name);};
        auto log_partial_path_of = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.partial_path_of,node2,node1,source_name);};
        auto log_conflicts_with  = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.conflicts_with,node2,node1,source_name);};
        auto log_resolved_by     = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.resolved_by,node2,node1,source_name);};
        auto log_terminal        = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.terminal,node2,node1,source_name);};

        perform_conflict_analysis(tree, ottid_to_node, mrca_of_pair,
                                  summaryTree, constSummaryOttIdToNode, mrca_of_pair,
                                  log_supported_by,
                                  log_partial_path_of,
                                  log_conflicts_with,
                                  log_resolved_by,
                                  log_terminal);
    }
    {
        auto nothing    = [](const node_t*, const node_t*) {};
        auto log_resolved_by     = [&source_name,&s](const node_t* node2, const node_t* node1) {add_element(s.resolves,node1,node2,source_name);};

        perform_conflict_analysis(summaryTree, constSummaryOttIdToNode, mrca_of_pair,
                                  tree, ottid_to_node, mrca_of_pair,
                                  nothing,
                                  nothing,
                                  nothing,
                                  log_resolved_by,
                                  nothing);
    }
}

void mapNextTree2(const Tree_t& summaryTree,
                  const map<OttId, const Tree_t::node_type*>& constSummaryOttIdToNode,
                  const Tree_t & tree,
                  stats& s) {
    typedef Tree_t::node_type node_type;
    std::function<const node_type*(const node_type*,const node_type*)> mrca_of_pair = [](const node_type* n1, const node_type* n2) {return mrca_from_depth(n1,n2);};
    string source_name = source_from_tree_name(summaryTree.get_name());
    auto ottid_to_node = get_ottid_to_const_node_map(tree);
    {
        auto log_supported_by    = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.supported_by,node2,node1,source_name);};
        auto log_partial_path_of = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.partial_path_of,node2,node1,source_name);};
        auto log_conflicts_with  = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.conflicts_with,node2,node1,source_name);};
        auto log_resolved_by     = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.resolved_by,node2,node1,source_name);};
        auto log_terminal        = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.terminal,node2,node1,source_name);};

        perform_conflict_analysis(tree, ottid_to_node, mrca_of_pair,
                                  summaryTree, constSummaryOttIdToNode, mrca_of_pair,
                                  log_supported_by,
                                  log_partial_path_of,
                                  log_conflicts_with,
                                  log_resolved_by,
                                  log_terminal);
    }
    {
        auto nothing    = [](const node_t*, const node_t*) {};
        auto log_resolved_by     = [&source_name,&s](const node_t* node1, const node_t* node2) {add_element(s.resolves,node1,node2,source_name);};
        perform_conflict_analysis(summaryTree, constSummaryOttIdToNode, mrca_of_pair,
                                  tree, ottid_to_node, mrca_of_pair,
                                  nothing,
                                  nothing,
                                  nothing,
                                  log_resolved_by,
                                  nothing);
    }
}

void mapNextTree(const Tree_t& summaryTree,
                 const map<OttId, const Tree_t::node_type*>& constSummaryOttIdToNode,
                 const Tree_t & tree,
                 stats& s, //isTaxoComp is third param
                 bool sw) {
    if (not sw) {
        mapNextTree1(summaryTree, constSummaryOttIdToNode, tree, s);
    } else {
        mapNextTree2(summaryTree, constSummaryOttIdToNode, tree, s);
    }
}

void show_header(std::ostream& o) {
    o << "supported_by" << "\t" << "partial_path_of" << "\t" << "terminal" << "\t" << "conflicts_with" << "\t" << "resolves" << "\t" << "resolved_by" << "\tname\n";
}

void show_stats(std::ostream& o, const stats& s, const string& name) {
    o << s.supported_by.size() << "\t" << s.partial_path_of.size() << "\t" << s.terminal.size() << "\t" << s.conflicts_with.size() << "\t" << s.resolves.size() << "\t" << s.resolved_by.size() << "\t" << name << "\n";
}

void show_group(std::ostream& o, const set<pair<string,string>>& group, const string& relation, const string& name) {
    for(auto& x: group) {
        o << "relation = " << relation << "   tree = " << x.first << "   node = " << x.second << "   group = " << name << "\n";
    }
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
        variables_map args = parse_cmd_line(argc, argv);
        vector<string> fir_astral = args["first"].as<vector<string> >();
        vector<string> sec_astral = args["second"].as<vector<string> >();
        ParsingRules rules;
        // 1. Load and process summary tree.
        auto fir_trees = get_trees<Tree_t>(fir_astral, rules);
        auto sec_trees = get_trees<Tree_t>(sec_astral, rules);
        // auto summaryOttIdToNode = get_ottid_to_node_map(*summaryTree);
        // auto constSummaryOttIdToNode = get_ottid_to_const_node_map(*summaryTree);
        // auto monotypic_nodes = suppress_and_record_monotypic(*summaryTree);
        // compute_depth(*summaryTree);
        // stats global;
        // // 2. Load and process input trees.
        // if (not names) {
        //     show_header(std::cout);
        // }
        // for(const auto& filename: inputs) {
        //     auto tree = get_tree<Tree_t>(filename);
        //     compute_depth(*tree);
        //     compute_summary_leaves(*tree, summaryOttIdToNode);
        //     string source_name = source_from_tree_name(tree->get_name());
        // }
    } catch (std::exception& e) {
        std::cerr << "otc-contrast-astral-runs: Error! " << e.what() << std::endl;
        exit(1);
    }
}
