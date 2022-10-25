#include <algorithm>
#include <set>
#include <list>
#include <iterator>
#include <numeric>

#include "solver/tree.h"
#include "solver/rsplit.h"
#include "solver/component.h"
#include "solver/solution.h"
#include "solver/build.h"
#include "solver/oracle.h"
#include "solver/combine.h"

#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
#include "otc/tree_iter.h"
#include "otc/induced_tree.h"
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>
#include <optional>

#include "otc/node_naming.h"

using namespace otc;
namespace fs = boost::filesystem;

using std::vector;
using std::unique_ptr;
using std::set;
using std::pair;
using std::list;
using std::map;
using std::string;
using std::optional;
using std::shared_ptr;

using namespace otc;

template <typename Tree_Out_t, typename Tree_In_t>
unique_ptr<Tree_Out_t> copy_tree(const Tree_In_t& tree)
{
    std::unique_ptr<Tree_Out_t> new_tree(new Tree_Out_t());

    // 1. Construct duplicate nodes for the new tree, recording correspondence
    std::unordered_map<const_node_type<Tree_In_t>*, non_const_node_type<Tree_Out_t>*> to_new_tree;
    for(auto nd: iter_post(tree))
    {
        auto nd2 = new_tree->create_node(nullptr);
        to_new_tree[nd] = nd2;

        if (nd->has_ott_id())
            nd2->set_ott_id(nd->get_ott_id());

        if (nd->get_name().size())
            nd2->set_name(nd->get_name());
    }

    // 2. Link corresponding nodes to their corresponding parents
    for(auto nd: iter_post(tree))
    {
        if (auto p = nd->get_parent())
        {
            auto nd2 = to_new_tree.at(nd);
            auto p2 = to_new_tree.at(p);
            p2->add_child(nd2);
        }
    }
    // 3. Set the root of the new tree to node corresponding to the MRCA
    new_tree->_set_root( to_new_tree.at(tree.get_root()) );
    return new_tree;
}

namespace po = boost::program_options;
using po::variables_map;

variables_map parse_cmd_line(int argc,char* argv[]) { 
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("subproblem", value<vector<string>>()->composing(),"File containing ordered subproblem trees.")
        ;

    options_description output("Standard options");
    output.add_options()
        ("incertae-sedis,I", value<string>(), "File containing Incertae sedis ids")
        ("root-name,n", value<string>(), "Rename the root to this name")
        ("no-higher-tips,l", "Tips may be internal nodes on the taxonomy.")
        ("prune-unrecognized,p","Prune unrecognized tips");

    options_description strategies("Solver strategies");
    strategies.add_options()
        ("batching",value<bool>()->default_value(false), "Make unresolved taxonomy from input tips.")
        ("oracle", value<bool>()->default_value(false), "Predict conflicting splits before BUILD.")
        ("incremental", value<bool>()->default_value(true),"Reuse work from previous BUILD.")
        ("rollback",value<bool>()->default_value(true), "Record rollback info in BUILDINC.")
        ;

    options_description other("Other options");
    other.add_options()
        ("synthesize-taxonomy,T","Make unresolved taxonomy from input tips.")
        ("allow-no-ids,a", "Allow problems w/o OTT ids")
        ("standardize,S", "Write out a standardized subproblem and exit.")
        ("input-deg-dist", value<string>(), "Write input trees degree distribution to filepath.")
        ("output-deg-dist", value<string>(), "Write output trees degree distribution to filepath.")
        ("time,m", "Report time taken to standard error.")
         ;

    options_description visible;
    visible.add(output).add(strategies).add(other).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("subproblem", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-solve-subproblem <trees-file1> [<trees-file2> ... ] [OPTIONS]\n"
                                                    "Takes a series of tree files.\n"
                                                    "Files are concatenated and the combined list treated as a single subproblem.\n"
                                                    "Trees should occur in order of priority, with the taxonomy last.",
                                                    visible, invisible, p);
    return vm;
}

struct Solution;

// A "partial" solution only has implied_splits + non_implied_splits
// A "full" solution also has components.

template <typename T>
bool is_subset(const std::set<T>& set_two, const std::set<T>& set_one) {
    return std::includes(set_one.begin(), set_one.end(), set_two.begin(), set_two.end());
}

optional<OttId> find_ancestor_id(const Tree_t::node_type* nd) {
    // Don't call this on the root node: it will abort.
    assert(nd->get_parent());
    while (auto p = nd->get_parent()) {
        if (p->has_ott_id()) {
            return p->get_ott_id();
        }
        nd = p;
    }
    return {};
}

bool is_ancestral_to(const Tree_t::node_type* anc, const Tree_t::node_type* n1) {
    if (depth(n1) < depth(anc)) {
        return false;
    }
    while(depth(n1) > depth(anc)) {
        assert(n1->get_parent());
        n1 = n1->get_parent();
    }
    assert(depth(n1) == depth(anc));
    return (n1 == anc);
}

map<const Tree_t::node_type*,const Tree_t::node_type*> check_placement(const Tree_t& summary,
                                                                       const Tree_t& taxonomy) {
    for(auto nd: iter_post(summary)) {
        if (nd->get_parent() and nd->get_name().size() and not nd->has_ott_id()) {
            LOG(WARNING)<<"Named taxonomy node has no OTTID!  Not checking for incertae sedis placement.";
            return {};
        }
    }
    map<const Tree_t::node_type*,const Tree_t::node_type*> placements;
    auto node_from_id = get_ottid_to_const_node_map(taxonomy);
    for(auto nd: iter_post(summary)) {
        if (nd->get_parent() and nd->has_ott_id()) {
            auto id = nd->get_ott_id();
            auto anc_id = find_ancestor_id(nd);
            if (not anc_id) { // ancestor is the root
                continue;
            }
            auto tax_nd = node_from_id.at(id);
            auto tax_anc = node_from_id.at(*anc_id);
            if (not is_ancestral_to(tax_anc, tax_nd)) {
                placements[tax_nd] = tax_anc;
            }
        }
    }
    return placements;
}


bool conflicting(const vector<int>& all_leaves_indices, const vector<ConstRSplit>& splits)
{
    return not BUILD_check(all_leaves_indices, splits);
}

bool conflicts_with(const vector<int>& all_leaves_indices, vector<ConstRSplit> splits1, const vector<ConstRSplit>& splits2)
{
    splits1.insert(splits1.end(), splits2.begin(), splits2.end());
    return conflicting(all_leaves_indices, splits1);
}

// Trying to find a conflicting set of splits where if you remove one split, then there is no longer a conflict.
// The first set of splits that conflicts?  Thinned from the back end to remove things that do not contribute to the conflict?
// Alternatively:
// Alternatively: the smallest set of splits that conflicts?

template<typename T>
std::vector<T> concat(const std::vector<T>& v1, const std::vector<T>& v2)
{
    auto v3 = v1;
    v3.insert(v3.end(), v2.begin(), v2.end());
    return v3;
}

void simplify(RSplitObj& s1, RSplitObj& s2)
{
    std::set<OttId> include1;
    std::set<OttId> exclude1;
    for(auto& i: s1.in)
        include1.insert(i);
    for(auto& o: s1.out)
        exclude1.insert(o);
    std::set<OttId> taxa1 = set_union_as_set(include1, exclude1);

    std::set<OttId> include2;
    std::set<OttId> exclude2;
    for(auto& i: s2.in)
        include2.insert(i);
    for(auto& o: s2.out)
        exclude2.insert(o);
    std::set<OttId> taxa2 = set_union_as_set(include2, exclude2);

    auto common_taxa = set_intersection_as_set(taxa1, taxa2);

    include1 = set_intersection_as_set(include1, common_taxa);
    exclude1 = set_intersection_as_set(exclude1, common_taxa);
    s1.in.clear();
    for(auto& i: include1)
        s1.in.push_back(i);
    s1.out.clear();
    for(auto& e: exclude1)
        s1.out.push_back(e);

    include2 = set_intersection_as_set(include2, common_taxa);
    exclude2 = set_intersection_as_set(exclude2, common_taxa);
    s2.in.clear();
    for(auto& i: include2)
        s2.in.push_back(i);
    s2.out.clear();
    for(auto& e: exclude2)
        s2.out.push_back(e);
}

// This is for debugging.
std::vector<ConstRSplit> find_minimal_conflict_set(const vector<int>& all_leaves_indices, const vector<ConstRSplit>& splits1, const vector<ConstRSplit>& splits2)
{
    assert(not conflicting(all_leaves_indices, splits1));
    assert(not conflicting(all_leaves_indices, splits2));
    assert(conflicts_with(all_leaves_indices, splits1, splits2));

    if (splits1.size() == 1)
    {
        return {splits1[0]};
    }
    int n_half = splits1.size()/2;

    // 1. If the first half conflicts, we can drop the last half.
    vector<ConstRSplit> splits1a;
    for(int i=0;i<n_half;i++)
        splits1a.push_back(splits1[i]);

    if (conflicts_with(all_leaves_indices, splits1a, splits2))
        return find_minimal_conflict_set(all_leaves_indices, splits1a, splits2);

    // 2. If the second half conflicts, we can drop the first half
    vector<ConstRSplit> splits1b;
    for(int i=n_half;i<splits1.size();i++)
        splits1b.push_back(splits1[i]);

    if (conflicts_with(all_leaves_indices, splits1b, splits2))
        return find_minimal_conflict_set(all_leaves_indices, splits1b, splits2);

    // 3. Find a minimal subset of splits1b to conflict with splits1a + splits2
    auto splits1b_conflicting = find_minimal_conflict_set(all_leaves_indices, splits1b, concat(splits1a, splits2));

    // 4. Find a minimal subset of splits1a to conflict with splits1b_conflicting + splits2;
    auto splits1a_conflicting = find_minimal_conflict_set(all_leaves_indices, splits1a, concat(splits1b_conflicting, splits2));

    return concat(splits1a_conflicting, splits1b_conflicting);
}

/// Create an unresolved taxonomy out of all the input trees.
unique_ptr<Tree_t> make_unresolved_tree(const vector<unique_ptr<Tree_t>>& trees, bool use_ids) {
    std::unique_ptr<Tree_t> retTree(new Tree_t());
    retTree->create_root();
    if (use_ids) {
        map<OttId,string> names;
        for(const auto& tree: trees) {
            for(auto nd: iter_pre(*tree)) {
                if (nd->is_tip()) {
                    OttId id = nd->get_ott_id();
                    auto it = names.find(id);
                    if (it == names.end()) {
                        names[id] = nd->get_name();
                    }
                }
            }
        }
        for(const auto& n: names) {
            auto node = retTree->create_child(retTree->get_root());
            node->set_ott_id(n.first);
            node->set_name(n.second);
        }
        clear_and_fill_des_ids(*retTree);
    } else {
        set<string> names;
        for(const auto& tree: trees) {
            for(auto nd: iter_pre(*tree)) {
                if (nd->is_tip()) {
                    names.insert(nd->get_name());
                }
            }
        }
        for(const auto& n: names) {
            auto node = retTree->create_child(retTree->get_root());
            node->set_name(n);
        }
    }
    return retTree;
}

string node_name_is(const node_t* nd, const OttIdSet & incertae_sedis) {
    std::ostringstream msg;
    if (incertae_sedis.count(nd->get_ott_id())) {
        msg << "?";
    }
    msg <<nd->get_name();
    if (nd->has_ott_id()) {
        msg << " [" << nd->get_ott_id() << "]";
    }
    return msg.str();
}


inline vector<const node_t *> vec_ptr_to_anc(const node_t * des, const node_t * anc) {
    vector<const node_t *> ret;
    while(des != anc) {
        assert(depth(des) > depth(anc));
        ret.push_back(des);
        des = des->get_parent();
    }
    ret.push_back(anc);
    return ret;
}

void standardize(Tree_t& t)
{
    std::unordered_map<const Tree_t::node_type*, OttId> smallest_child;
    calculate_smallest_child_map<Tree_t>(t, smallest_child);
    sort_by_smallest_child_map(t, smallest_child);
}

int main(int argc, char *argv[])
{
    std::cout<<std::boolalpha;
    std::cerr<<std::boolalpha;
    try {
        // 1. Parse command line arguments
        variables_map args = parse_cmd_line(argc,argv);
        ParsingRules rules;
        rules.set_ott_ids = not (bool)args.count("allow-no-ids");
        rules.prune_unrecognized_input_tips = (bool)args.count("prune-unrecognized");
        bool synthesize_taxonomy = (bool)args.count("synthesize-taxonomy");
        bool cladeTips = not (bool)args.count("no-higher-tips");
        bool verbose = (bool)args.count("verbose");
        bool writeStandardized = (bool)args.count("standardize");
        if (writeStandardized) {
            rules.set_ott_ids = false;
        }
        bool setRootName = (bool)args.count("root-name");
        vector<string> filenames = args["subproblem"].as<vector<string>>();
        // 2. Load trees from subproblem file(s)
        if (filenames.empty()) {
            throw OTCError("No subproblem provided!");
        }
        vector<unique_ptr<Tree_t>> trees = get_trees<Tree_t>(filenames, rules);
        if (trees.empty()) {
            throw OTCError("No trees loaded!");
        }
        if (args.count("input-deg-dist")) {
            auto filename = args["input-deg-dist"].as<string>();
            std::ofstream file(filename); 
            if (not file) {
                throw OTCError() << "Cannot open input-deg-dist file '" << fs::absolute(filename) << "'";
            }
            for(const auto & tree: trees) {
                writeDegDist(file, nullptr, *tree);
            }
        }
        //2.5 Load Incertae Sedis info
        OttIdSet incertae_sedis;
        if (args.count("incertae-sedis")) {
            auto filename = args["incertae-sedis"].as<string>();
            std::ifstream file(filename);
            if (not file) {
                throw OTCError() << "Cannot open incertae sedis file '" << fs::absolute(filename) << "'";
            }
            while (file) {
                OttId i;
                file >> i;
                incertae_sedis.insert(i);
            }
        }
        // 3. Make a fake taxonomy if asked
        if (synthesize_taxonomy) {
            trees.push_back(make_unresolved_tree(trees, rules.set_ott_ids));
            LOG(DEBUG) << "taxonomy = " << newick(*trees.back()) << "\n";
        }
        // 4. Add fake Ott Ids to tips and compute des_ids (if asked)
        if (not rules.set_ott_ids) {
            auto name_to_id = create_ids_from_names(*trees.back());
            for(auto& tree: trees) {
                set_ids_from_names_and_refresh(*tree, name_to_id);
            }
        }
        // 5. Write out subproblem with newly minted ottids (if asked)
        if (writeStandardized) {
            for(const auto& tree: trees) {
                relabel_nodes_with_ott_id(*tree);
                std::cout << newick(*tree) << "\n";
            }
            return 0;
        }
        // 6. Check if trees are mapping to non-terminal taxa, and either fix the situation or die.
        for (int i = 0; i < trees.size() - 1; i++) {
            if (cladeTips) {
                expand_ott_internals_which_are_leaves(*trees[i], *trees.back());
            } else {
                require_tips_to_be_mapped_to_terminal_taxa(*trees[i], *trees.back());
            }
        }
        // 6.5 Make a copy of the taxonomy so that "combine" doesn't modify it.
        auto taxonomy = copy_tree<Tree_t>(*trees.back());
        compute_depth(*taxonomy);

        // 7. Perform the synthesis
        auto tree = combine(trees, incertae_sedis, args);
        // 8. Set the root name (if asked)
        // FIXME: This could be avoided if the taxonomy tree in the subproblem always had a name for the root node.
        if (setRootName) {
            tree->get_root()->set_name(args["root-name"].as<string>());
        }
        // 9. Write out the summary tree.
        standardize(*tree);
        write_tree_as_newick(std::cout, *tree);
        std::cout << "\n";

        if (args.count("output-deg-dist")) {
            auto filename = args["output-deg-dist"].as<string>();
            std::ofstream file(filename); 
            if (not file) {
                throw OTCError() << "Cannot open output-deg-dist file '" << fs::absolute(filename) << "'";
            }
            writeDegDist(file, nullptr, *tree);
        }

        // 10. Find placements
        auto placements = check_placement(*tree, *taxonomy);
        for(auto& [placed, parent]: placements)
        {
            auto mrca = mrca_from_depth(placed, parent);
            vector<const node_t*> placement_path = vec_ptr_to_anc(parent, mrca);
            vector<const node_t*> is_path = vec_ptr_to_anc(placed, mrca);
            std::ostringstream msg;
            for(int i=0;i<is_path.size();i++) {
                msg <<node_name_is(is_path[i], incertae_sedis);
                if (i != is_path.size()-1) {
                    msg << " <- ";
                }
            }
            msg << " placed under ";
            for(int i=0;i<placement_path.size();i++) {
                msg << node_name_is(placement_path[i], incertae_sedis);
                if (i != placement_path.size()-1) {
                    msg << " <- ";
                }
            }
            LOG(INFO) << msg.str();
        }
        return 0;
    } catch (std::exception& e) {
        std::cerr << "otc-solve-subproblem: Error! " << e.what() << std::endl;
        exit(1);
    }
}
