#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"

#include <boost/range/adaptor/reversed.hpp>

INITIALIZE_EASYLOGGINGPP

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unique_ptr;

using boost::spirit::qi::symbols;
using boost::string_ref;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;

variables_map parse_cmd_line(int argc,char* argv[]) 
{ 
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("tree", value<string>(),"Newick file for tree")
        ;

    options_description tree("Tree options");
    tree.add_options()
        ("root,r", value<long>(), "OTT id of root node of subtree to keep")
        ("prune,p", value<std::vector<long> >()->composing(),"OTT ids of taxa to prune")
        ("slice,s", value<std::vector<long> >()->composing(),"OTT ids of root and taxa to prune")
        ;

    options_description output("Output options");
    output.add_options()
        ("high-degree-nodes", value<long>(), "Show the top <arg> high-degree nodes.")
        ("degree-of",value<long>(), "Show the degree of node <arg>")
        ("children-of",value<long>(), "List the children of node <arg>")
        ("parent-of",value<long>(), "List the parent of node <arg>")
        ("count-nodes","Show the number of leaves")
        ("count-leaves","Show the number of leaves")
        ;

    options_description visible;
    visible.add(tree).add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("tree", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: tree-tool <newick-tree> [OPTIONS]\n"
                                                    "Do various operations on newick trees, such as extracting subtrees by OTT id.",
                                                    visible, invisible, p);
    return vm;
}

long n_nodes(const Tree_t& T) {
#pragma clang diagnostic ignored  "-Wunused-variable"
#pragma GCC diagnostic ignored  "-Wunused-variable"
    long count = 0;
    for(auto nd: iter_post_const(T)){
        count++;
    }
    return count;
}

long n_leaves(const Tree_t& T) {
#pragma clang diagnostic ignored  "-Wunused-variable"
#pragma GCC diagnostic ignored  "-Wunused-variable"
    long count = 0;
    for(auto nd: iter_leaf_const(T)){
        count++;
    }
    return count;
}

unique_ptr<Tree_t> get_tree(const string& filename)
{
    vector<unique_ptr<Tree_t>> trees;
    std::function<bool(unique_ptr<Tree_t>)> a = [&](unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;};
    ParsingRules rules;
    rules.requireOttIds = false;
    otc::processTrees(filename,rules,a);//[&](unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;});
    return std::move(trees[0]);
}

Tree_t::node_type* find_node_by_ott_id(Tree_t& tree, long root_ott_id)
{
    for(auto nd: iter_pre(tree))
        if (nd->hasOttId() and nd->getOttId() == root_ott_id)
            return nd;
    
    throw OTCError()<<"Can't find node with id "<<root_ott_id<<" in tree '"<<tree.getName()<<"'";
}

Tree_t::node_type* find_node_by_name(Tree_t& tree, const string& name)
{
    for(auto nd: iter_pre(tree))
        if (nd->getName().size() and nd->getName() == name)
            return nd;
    
    throw OTCError()<<"Can't find node with name '"<<name<<"' in tree '"<<tree.getName()<<"'";
}

unique_ptr<Tree_t> truncate_to_subtree_by_ott_id(unique_ptr<Tree_t> tree, long root_ott_id)
{
    auto root = find_node_by_ott_id(*tree, root_ott_id);
    root->detachThisNode();
    unique_ptr<Tree_t> tree2 (new Tree_t);
    tree2->_setRoot(root);
    return tree2;
}

unique_ptr<Tree_t> slice_tree(unique_ptr<Tree_t> tree,
                              long root_ott_id,
                              const std::vector<long> & tips) {
    for (auto t_ott_id : tips) {
        auto tn = find_node_by_ott_id(*tree, t_ott_id);
        if (tn && !tn->isTip()) {
            const auto children = all_children(tn);
            for (auto c : children) {
                c->detachThisNode();
            }
        }
    }
    auto root = find_node_by_ott_id(*tree, root_ott_id);
    root->detachThisNode();
    unique_ptr<Tree_t> tree2 (new Tree_t);
    tree2->_setRoot(root);
    return tree2;
}


void show_high_degree_nodes(const Tree_t& tree, int n)
{
    std::multimap<long,std::string> nodes;
    for(auto nd: iter_pre_const(tree))
    {
        auto outdegree = nd->getOutDegree();
        if (outdegree > 1)
            nodes.insert({outdegree,nd->getName()});
    }

    for(const auto& x: boost::adaptors::reverse(nodes))
    {
        --n;
        if (n <= 0) return;
        std::cout<<x.first<<"\t"<<x.second<<std::endl;
    } 
}

int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);

    try
    {
        variables_map args = parse_cmd_line(argc,argv);

        if (not args.count("tree"))
            throw OTCError()<<"Please specify the tree to operate on!";
        
        auto tree = get_tree(args["tree"].as<string>());

        if (args.count("root"))
        {
            long root = args["root"].as<long>();
            tree = truncate_to_subtree_by_ott_id(std::move(tree), root);
        }
        if (args.count("slice"))
        {
            std::vector<long> slice = args["slice"].as<std::vector<long> >();
            if (slice.empty()) {
                throw OTCError() << "Expecting a root ID followed by a OTT Ids to slice from the tree";
            }
            long root = slice[0];
            slice.erase(slice.begin());
            tree = slice_tree(std::move(tree), root, slice);
        }
        
        if (args.count("high-degree-nodes"))
        {
            long n = args["high-degree-nodes"].as<long>();
            show_high_degree_nodes(*tree, n);
        }
        else if (args.count("degree-of"))
        {
            long n = args["degree-of"].as<long>();
            auto nd = find_node_by_ott_id(*tree, n);
            std::cout<<nd->getOutDegree()<<"\n";
        }
        else if (args.count("children-of"))
        {
            long n = args["children-of"].as<long>();
            auto nd = find_node_by_ott_id(*tree, n);
            for(auto c = nd->getFirstChild(); c; c = c->getNextSib())
                std::cout<<c->getName()<<"\n";
        }
        else if (args.count("parent-of"))
        {
            long n = args["parent-of"].as<long>();
            auto nd = find_node_by_ott_id(*tree, n);
            if (nd->getParent())
                std::cout<<nd->getParent()->getName()<<"\n";
            else
                std::cout<<"No parent: that node is the root.\n";
        }
        else if (args.count("count-nodes"))
        {
            std::cout<<n_nodes(*tree)<<std::endl;
        }
        else if (args.count("count-leaves"))
        {
            std::cout<<n_leaves(*tree)<<std::endl;
        }
        else {
            writeTreeAsNewick(std::cout, *tree);
            std::cout << std::endl;
        }
    }
    catch (std::exception& e)
    {
        cerr<<"otc-tree-tool: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}

// 1. Write a parser to read the lines faster
// 2. Avoid memory allocation -- by mmapping the taxonomy file?
// 3. Convert the flags into a bitmask
// 4. Should the Rank be a converted to an integer?
// 5. Can we assign OTT IDs to internal nodes of a tree while accounting for Incertae Sedis taxa?
// * What are the triplet-inference rules for the Incertae Sedis problem?

// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

