#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <regex>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "otc/node_naming.h"

#include <boost/range/adaptor/reversed.hpp>

namespace fs = std::filesystem;

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unique_ptr;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeSmallestChild, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;

variables_map parse_cmd_line(int argc,char* argv[]) { 
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("tree", value<string>(),"Newick file for tree")
        ;

    options_description tree("Tree options");
    tree.add_options()
        ("root,r", value<OttId>(), "OTT id of root node of subtree to keep")
        ("prune,p", value<std::vector<OttId> >()->composing(),"OTT ids of taxa to prune")
        ("slice,s", value<std::vector<OttId> >()->composing(),"OTT ids of root and taxa to prune")
        ;

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("taxonomy", value<string>(),"Directory name for the taxonomy")
        ;

    options_description output("Output options");
    output.add_options()
        ("high-degree-nodes", value<long>(), "Show the top <arg> high-degree nodes.")
        ("degree-of",value<OttId>(), "Show the degree of node <arg>")
        ("children-of",value<OttId>(), "List the children of node <arg>")
        ("parent-of",value<OttId>(), "List the parent of node <arg>")
        ("ancestors-of",value<OttId>(), "List the all ancestor (parent to root order) of node <arg>")
        ("count-nodes","Show the number of nodes")
        ("count-leaves","Show the number of leaves")
        ("show-leaves","Show the number of leaves")
        ("show-internal","Show the number of leaves")
        ("write-taxonomy",value<string>(),"Write as taxonomy in directory <arg>")
        ("lost-taxa-vs",value<string>(),"Taxonomy tree to compare for lost taxa.")
        ("standardize","Perform a rotation to a standard form.")
        ("indented-table","print number of leaves for each internal node")
        ;

    options_description visible;
    visible.add(tree).add(taxonomy).add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("tree", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: tree-tool <newick-tree> [OPTIONS]\n"
                                                    "Do various operations on newick trees, such as extracting subtrees by OTT id.",
                                                    visible, invisible, p);
    return vm;
}


void show_nodes(const Tree_t& T) {
    for(auto nd: iter_post_const(T)) {
        if (nd->get_name().size()) {
            std::cout << nd->get_name() << "\n";
        }
    }
}

std::size_t n_leaves(const Tree_t& T) {
#pragma clang diagnostic ignored  "-Wunused-variable"
#pragma GCC diagnostic ignored  "-Wunused-variable"
    std::size_t count = 0;
    for(auto nd: iter_leaf_const(T)){
        count++;
    }
    return count;
}

void show_leaf_nodes(const Tree_t& T) {
    for(auto nd: iter_leaf_const(T)) {
        if (nd->get_name().size()){
            std::cout << nd->get_name() << "\n";
        }
    }
}

void show_internal_nodes(const Tree_t& T) {
    for(auto nd: iter_post_const(T)) {
        if (nd->get_name().size() and not nd->is_tip()) {
            std::cout << nd->get_name() << "\n";
        }
    }
}

void indented_table_of_node_counts(std::ostream & out, const Tree_t & tree) {
    std::map<const Tree_t::node_type*, long> nd2numLeaves;
    std::map<const Tree_t::node_type*, long> nd2Indent;
    long count = 0;
    for(auto nd: iter_post_const(tree)){
        auto p = nd->get_parent();
        if (p) {
            if (nd->is_tip()) {
                nd2numLeaves[p] += 1;
            } else {
                nd2numLeaves[p] += nd2numLeaves[nd];
            }
        }
    }
    for(auto nd: iter_pre_const(tree)){
        if (nd->is_tip()) {
            continue;
        }
        auto p = nd->get_parent();
        if (p) {
            nd2Indent[nd] = 2 + nd2Indent[p];
        } else {
            nd2Indent[nd] = 0;
        }
        out << std::string(nd2Indent[nd], ' ');
        out << nd->get_name() << " : " << nd2numLeaves[nd] << '\n';
    }   
    
}

unique_ptr<Tree_t> get_tree(const string& filename) {
    vector<unique_ptr<Tree_t>> trees;
    std::function<bool(unique_ptr<Tree_t>)> a = [&](unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;};
    ParsingRules rules;
    rules.require_ott_ids = false;
    otc::process_trees(filename,rules,a);//[&](unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;});
    return std::move(trees[0]);
}

Tree_t::node_type* find_node_by_ott_id(Tree_t& tree, OttId root_ott_id, bool throw_if_not_found=true) {
    for(auto nd: iter_pre(tree)) {
        if (nd->has_ott_id() and nd->get_ott_id() == root_ott_id) {
            return nd;
        }
    }
    if (throw_if_not_found) {
        throw OTCError() << "Can't find node with id " << root_ott_id << " in tree '" << tree.get_name() << "'";
    }
    return nullptr;
}

Tree_t::node_type* find_node_by_name(Tree_t& tree, const string& name) {
    for(auto nd: iter_pre(tree)) {
        if (nd->get_name().size() and nd->get_name() == name) {
            return nd;
        }
    }
    throw OTCError() << "Can't find node with name '" << name << "' in tree '" << tree.get_name() << "'";
}

unique_ptr<Tree_t> truncate_to_subtree_by_ott_id(unique_ptr<Tree_t> tree, OttId root_ott_id) {
    auto root = find_node_by_ott_id(*tree, root_ott_id);
    root->detach_this_node();
    unique_ptr<Tree_t> tree2 (new Tree_t);
    tree2->_set_root(root);
    return tree2;
}

void prune_from_tree(Tree_t & tree, const std::vector<OttId> & tips) {
    OttIdSet tipset;
    for (auto t : tips) {
        tipset.insert(t);
    }
    std::set<Tree_t::node_type *> todel;
    for(auto nd: iter_pre(tree)) {
        if (nd->has_ott_id() and tipset.count(nd->get_ott_id()) > 0) {
            todel.insert(nd);
        }
    }
    for (auto tdn : todel) {
        tdn->detach_this_node();
    }
}

void prune_tree(Tree_t & tree, const std::vector<OttId> & tips) {
    for (auto t_ott_id : tips) {
        auto tn = find_node_by_ott_id(tree, t_ott_id);
        if (tn && !tn->is_tip()) {
            const auto children = all_children(tn);
            for (auto c : children) {
                c->detach_this_node();
            }
        }
    }
}

unique_ptr<Tree_t> slice_tree(unique_ptr<Tree_t> tree,
                              OttId root_ott_id,
                              const std::vector<OttId> & tips) {
    prune_tree(*tree, tips);
    auto root = find_node_by_ott_id(*tree, root_ott_id);
    root->detach_this_node();
    unique_ptr<Tree_t> tree2 (new Tree_t);
    tree2->_set_root(root);
    return tree2;
}


void show_high_degree_nodes(const Tree_t& tree, int n) {
    std::multimap<OttId,std::string> nodes;
    for(auto nd: iter_pre_const(tree)) {
        auto outdegree = nd->get_out_degree();
        if (outdegree > 1) {
            nodes.insert({outdegree,nd->get_name()});
        }
    }
    for(const auto& x: boost::adaptors::reverse(nodes)) {
        --n;
        if (n <= 0) {
            return;
        }
        std::cout << x.first << "\t" << x.second << std::endl;
    }
}

void create_file( const fs::path & ph, const std::string & contents ) {
    std::ofstream f( ph.string().c_str() );
    if (not f) {
        throw OTCError() << "Could not create empty file '" << ph.string() << "'";
    }
    if (not contents.empty()) {
        f << contents;
    }
}

std::string remove_ott_suffix(std::string name) {
    static std::regex ott("(.*)[_ ]ott.*");
    std::smatch matches;
    if (std::regex_match(name,matches,ott)) {
        name = matches[1];
    }
    return name;
}

void writeTreeAsTaxonomy(const string& dirname, const Tree_t& tree) {
    fs::path new_dir = dirname;
    if (! fs::exists(new_dir)) {
        fs::create_directories(new_dir);
    }
    // Create empty files
    for(const auto& name: {"conflicts.tsv", "deprecated.tsv", "log.tsv", "otu_differences.tsv", "weaklog.csv"}) {
        create_file(new_dir / name, "");
    }
    // Write the about.json file.
    create_file(new_dir/"about.json",string("{\"inputs\":[\"") + tree.get_name() + "\"]}");
    // Write the synonyms.tsv file.
    create_file(new_dir/"synonyms.tsv","name\t|\tuid\t|\ttype\t|\tuniqname\t|\t\n");
    // Write the version file.
    create_file(new_dir/"version.txt","0.0");
    // Write the new taxonomy file.
    {
        std::ofstream tf ((new_dir/"taxonomy.tsv").string());
        tf << "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t" << std::endl;
        string sep = "\t|\t";
        for(auto nd: iter_pre_const(tree)) {
            tf << nd->get_ott_id();
            tf << sep;
            if (nd->get_parent()) {
                tf << nd->get_parent()->get_ott_id();
            }
            tf << sep;
            tf << remove_ott_suffix(nd->get_name());
            tf << sep;
            tf << "no rank";
            tf << sep;
            tf << "tree:0";
            //      tf << "tree:" << tree.get_name();
            tf << sep;
            tf << sep;
            tf << sep;
            tf << "\n";
        }
        tf.close();
    }
    // Write the new forwards file.
    {
        std::ofstream ff((new_dir/"forwards.tsv").string());
        ff << "id\treplacement\n";
        ff.close();
    }
}

// BDR - 7/17/2017
// There is similar code in tools/taxonomy-parser.cpp.  ALthough not that similar.
// Maybe we should merge the two.  This would require adding a --tax-root argument, and changing
//  the load_taxonomy(args) function in taxonomy.h to understand it.
void show_lost_taxa(const Tree_t& tree, const string& tax_tre_filename, const string& tax_tax_filename) {
    Taxonomy taxonomy{tax_tax_filename};
    auto tax_tree = first_newick_tree_from_file<Tree_t>(tax_tre_filename);
    auto all_taxa = get_all_ott_ids(*tax_tree);
    auto unbroken_taxa = get_all_ott_ids(tree);
    auto broken_taxa = set_difference_as_set(all_taxa, unbroken_taxa);
    vector<const TaxonomyRecord*> records;
    for(auto& id: broken_taxa) {
        records.push_back( &taxonomy.record_from_id(id) );
    }
    std::sort(records.begin(), records.end(), [](const auto& a, const auto& b) {return a->depth < b->depth;});
    for(const auto& rec: records) {
        std::cout << "depth=" << rec->depth << "   id=" << rec->id << "   uniqname='" << rec->uniqname << "'\n";
    }
}


void standardize(Tree_t& tree)
{
    calculate_smallest_child(tree);
    sort_by_smallest_child(tree);
}

int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);
    try {
        variables_map args = parse_cmd_line(argc,argv);
        if (not args.count("tree")) {
            throw OTCError() << "Please specify the tree to operate on!";
        }
        auto tree = get_tree(args["tree"].as<string>());
        if (args.count("root")) {
            OttId root = args["root"].as<OttId>();
            tree = truncate_to_subtree_by_ott_id(std::move(tree), root);
        }
        if (args.count("slice")) {
            std::vector<OttId> slice = args["slice"].as<std::vector<OttId> >();
            if (slice.empty()) {
                throw OTCError() << "Expecting a root ID followed by a OTT Ids to slice from the tree";
            }
            OttId root = slice[0];
            slice.erase(slice.begin());
            tree = slice_tree(std::move(tree), root, slice);
        }
        if (args.count("prune")) {
            std::vector<OttId> prune = args["prune"].as<std::vector<OttId> >();
            if (prune.empty()) {
                throw OTCError() << "OTT Ids to prune from the tree";
            }
            prune_from_tree(*tree, prune);
        }
        if (args.count("high-degree-nodes")) {
            long n = args["high-degree-nodes"].as<long>();
            show_high_degree_nodes(*tree, n);
        } else if (args.count("degree-of")) {
            OttId n = args["degree-of"].as<OttId>();
            auto nd = find_node_by_ott_id(*tree, n);
            std::cout << nd->get_out_degree() << "\n";
        } else if (args.count("children-of")) {
            OttId n = args["children-of"].as<OttId>();
            auto nd = find_node_by_ott_id(*tree, n);
            for(auto c = nd->get_first_child(); c; c = c->get_next_sib()) {
                std::cout << c->get_name() << "\n";
            }
        } else if (args.count("parent-of")) {
            OttId n = args["parent-of"].as<OttId>();
            auto nd = find_node_by_ott_id(*tree, n);
            if (nd->get_parent()) {
                std::cout << nd->get_parent()->get_name() << "\n";
            } else {
                std::cout << "No parent: that node is the root.\n";
            }
        } else if (args.count("ancestors-of")) {
            OttId n = args["ancestors-of"].as<OttId>();
            auto nd = find_node_by_ott_id(*tree, n);
            auto ancnd = nd->get_parent();
            if (ancnd) {
                while (ancnd) {
                    std::cout << ancnd->get_name() << "\n";
                    ancnd = ancnd->get_parent();
                }
            } else {
                std::cout << "No parent: that node is the root.\n";
            }
        } else if (args.count("count-nodes")) {
            std::cout << n_nodes(*tree) << std::endl;
        } else if (args.count("show-nodes")) {
            show_nodes(*tree);
        } else if (args.count("count-leaves")) {
            std::cout << n_leaves(*tree) << std::endl;
        } else if (args.count("show-leaves")) {
            show_leaf_nodes(*tree);
        } else if (args.count("show-internal")) {
            show_internal_nodes(*tree);
        } else if (args.count("indented-table")) {
            indented_table_of_node_counts(std::cout, *tree);
        } else if (args.count("write-taxonomy")) {
            string dirname = args["write-taxonomy"].as<string>();
            writeTreeAsTaxonomy(dirname, *tree);
        } else if (args.count("lost-taxa-vs")) {
            string tax_tre_filename = args["lost-taxa-vs"].as<string>();
            string tax_tax_filename = args["taxonomy"].as<string>();
            show_lost_taxa(*tree, tax_tre_filename, tax_tax_filename);
        } else {
            if (args.count("standardize"))
                standardize(*tree);
            write_tree_as_newick(std::cout, *tree);
            std::cout << std::endl;
        }
    } catch (std::exception& e) {
        cerr << "otc-tree-tool: Error! " << e.what() << std::endl;
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

