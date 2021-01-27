#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <regex>
#include <utility>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/newick.h"
#include "otc/tree_operations.h"
#include "otc/induced_tree.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"

#include <boost/range/adaptor/reversed.hpp>

#include <boost/filesystem/operations.hpp>

#include "otc/ws/prune.h"

namespace fs = boost::filesystem;

INITIALIZE_EASYLOGGINGPP

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unique_ptr;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

struct RTNodeDepth {
    int depth = 0; // depth = number of nodes to the root of the tree including the  endpoints (so depth of root = 1)
};

using Tree_t = RootedTree<RTNodeDepth, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;

variables_map parse_cmd_line(int argc,char* argv[]) { 
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("trees", value<vector<string>>()->composing(),"Filenames for newick trees")
        ;

    options_description output("Output options");
    output.add_options()
        ("out-dir",value<string>(),"Output directory for the newick files")
        ;

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("config,c",value<string>(),"Config file containing flags to filter")
        ("clean",value<string>(),"Comma-separated string of flags to filter")
        ("root,r", value<OttId>(), "OTT id of root node of subtree to keep")
        ("taxonomy", value<string>(),"Directory name for the taxonomy")
        ;


    options_description visible;
    visible.add(output).add(taxonomy).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("trees", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: prune-trees <newick-tree1>:<name1> <newick-tree2>:<name2> ... [OPTIONS]\n"
                                                    "Prune flagged taxa and remove unmapped tips, writing resulting Newick files to out-dir.",
                                                    visible, invisible, p);
    return vm;
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

std::string remove_ott_suffix(std::string name) {
    static std::regex ott("(.*)[_ ]ott.*");
    std::smatch matches;
    if (std::regex_match(name,matches,ott)) {
        name = matches[1];
    }
    return name;
}

void write_tree(const Tree_t& tree, const fs::path& out_path)
{
    std::ofstream out_file( out_path.string() );
    if (not out_file)
    {
        throw OTCError() << "Could not create empty file '" << out_path.string() << "'";
    }
    write_tree_as_newick(out_file, tree);
}


// OK, so we need to make an induced tree from a taxonomy.

// Remove leaves (and their monotypic ancestors) in the query tree
// that are the ancestors of other leaves in the query tree.

void prune_ancestral_leaves(Tree_t& query_tree, const Tree_t& taxonomy_tree, const std::unordered_map<OttId,const Tree_t::node_type*>& tax_id_to_node)
{
    // 1. Get an induced taxonomy tree starting from the leaves of the query tree
    auto query_id_to_node = get_ottid_to_const_node_map(query_tree);
    auto induced_leaves = get_induced_leaves<const Tree_t,const Tree_t>(taxonomy_tree, tax_id_to_node, query_tree, query_id_to_node);
    auto mrca = [](const Tree_t::node_type* n1, const Tree_t::node_type* n2) {return mrca_from_depth(n1,n2);};
    auto induced_taxonomy = get_induced_tree<Tree_t>(induced_leaves, mrca);

    // 2. Find the query tree leaves that are not tips on the induced taxonomy
    auto ottid_to_induced_tax_node = get_ottid_to_const_node_map(*induced_taxonomy);

    vector<Tree_t::node_type*> nodes_to_prune;
    for(auto leaf: iter_leaf(query_tree))
    {
        auto ottid = leaf->get_ott_id();
        auto tax_node = ottid_to_induced_tax_node.at(ottid);
        if (not tax_node->is_tip()) {
            nodes_to_prune.push_back(leaf);
        }
    }

    // 3. Prune the leaves that we selected while walking the query tree
    for(auto node: nodes_to_prune)
        delete_tip_and_monotypic_ancestors(query_tree, node);
}


std::pair<string,string> split_on_last(const string& s, char c)
{
    auto pos = s.rfind(c);
    if (pos == string::npos)
        return {s,""};
    else
        return {s.substr(0,pos),s.substr(pos+1)};
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    try {
        variables_map args = parse_cmd_line(argc,argv);

        // Where does this load the taxonomy from?
        // Should I remove the --config argument?  Or should we pass the propinquity config here?
        LOG(INFO)<<"Loading taxonomy:";
        auto taxonomy = load_taxonomy(args);
        LOG(INFO)<<"  "<<taxonomy.get_version();

        LOG(INFO)<<"Making taxonomy tree:";
        auto taxonomy_tree = taxonomy.get_tree<Tree_t>([](auto&){return "";});
        compute_depth(*taxonomy_tree);
        auto tax_node_map = get_ottid_to_const_node_map(*taxonomy_tree);
        LOG(INFO)<<"done.";

        if (not args.count("trees")) throw OTCError() << "No trees given!";

        auto filenames = args["trees"].as<vector<string>>();

        if (not args.count("out-dir"))
            throw OTCError()<<"output directory not specified!  Use --out-dir=<directory>";
            
        auto out_dir = fs::path(args["out-dir"].as<string>());

        for(auto& filename_name: filenames)
        {
            // We should split this on ':', and use the last part as a name to write out.
            auto [in_filename, out_name] = split_on_last(filename_name, ':');
            if (out_name.empty())
                throw OTCError()<<"tree file '"<<filename_name<<"' should have the form <filename>:<name>";

            auto tree = get_tree(in_filename);
            compute_depth(*tree);

            int n_leaves1 = n_leaves(*tree);

            prune_unmapped_leaves(*tree, taxonomy);

            // Clear the existing names to generate ott_XXXX names from the ids,
            // since the ids have just been forwarded by prune_unmapped_leaves( ).
            for(auto leaf: iter_leaf(*tree))
            {
                assert(leaf->has_ott_id());
                leaf->set_name("");
            }

            int n_leaves2 = n_leaves(*tree);

            prune_duplicate_ottids(*tree);

            int n_leaves3 = n_leaves(*tree);

            prune_ancestral_leaves(*tree, *taxonomy_tree, tax_node_map);

            int n_leaves4 = n_leaves(*tree);

            // Uh... what tree are we supposed to write here?
            // write_tree(*tree, out_dir / (out_name + "-taxonomy.tre"));

            // Write out the pruned tree
            write_tree(*tree, out_dir / (out_name + ".tre"));

            LOG(INFO)<<"Pruning tree '"<<out_name<<"' from "<<n_leaves1<<" to "<<n_leaves4<<" leaves.";
            LOG(INFO)<<"   no OTT ID: "<<n_leaves1 - n_leaves2<<"   duplicate OTT id: "<<n_leaves2 - n_leaves3<<"   leaf ancestor of other leaf: "<<n_leaves3 - n_leaves4;
        }
    } catch (std::exception& e) {
        cerr << "otc-prune-trees: Error! " << e.what() << std::endl;
        exit(1);
    }
}
