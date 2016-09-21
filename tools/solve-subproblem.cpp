#include <algorithm>
#include <set>
#include <list>
#include <iterator>

#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
#include "otc/tree_iter.h"
#include <fstream>

using namespace otc;
using std::vector;
using std::unique_ptr;
using std::set;
using std::list;
using std::map;
using std::string;
using namespace otc;

typedef TreeMappedWithSplits Tree_t;

/// Create a SORTED vector from a set
template <typename T>
vector<T> set_to_vector(const set<T>& s) {
    vector<T> v;
    v.reserve(s.size());
    std::copy(s.begin(), s.end(), std::back_inserter(v));
    return v;
}

struct RSplit {
    vector<int> in;
    vector<int> out;
    vector<int> all;
    RSplit() = default;
    RSplit(const set<int>& i, const set<int>& a)
        {
            in  = set_to_vector(i);
            all = set_to_vector(a);
            set_difference(begin(all), end(all), begin(in), end(in), std::inserter(out, out.end()));
            assert(in.size() + out.size() == all.size());
        }
};

std::ostream& operator<<(std::ostream& o, const RSplit& s);
int merge_components(int c1, int c2, vector<int>& component, vector<list<int>>& elements);
bool empty_intersection(const set<int>& xs, const vector<int>& ys);
unique_ptr<Tree_t> BUILD(const std::vector<int>& tips, const std::vector<const RSplit*>& splits);
unique_ptr<Tree_t> BUILD(const vector<int>& tips, const vector<RSplit>& splits);
void add_names(unique_ptr<Tree_t>& tree, const unique_ptr<Tree_t>& taxonomy);
set<int> remap_ids(const set<long>& s1, const map<long,int>& id_map);
unique_ptr<Tree_t> combine(const vector<unique_ptr<Tree_t> >& trees, const set<long>&, bool verbose);
unique_ptr<Tree_t> make_unresolved_tree(const vector<unique_ptr<Tree_t>>& trees, bool use_ids);

namespace po = boost::program_options;
using po::variables_map;

variables_map parse_cmd_line(int argc,char* argv[]) 
{ 
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("subproblem", value<vector<string>>()->composing(),"File containing ordered subproblem trees.")
        ;

    options_description output("Some options");
    output.add_options()
        ("standardize,S", "Write out a standardized subproblem and exit.")
        ("synthesize-taxonomy,T","Synthesize an unresolved taxonomy from all mentioned tips.")
        ("no-higher-tips,l", "Tips may be internal nodes on the taxonomy.")
        ("root-name,n",value<string>(), "Rename the root to this name")
        ("require-ott-ids,o", "Require OTT ids")
        ("prune-unrecognized,p","Prune unrecognized tips")
	("incertae-sedis,I",value<string>(),"File containing Incertae sedis ids")
        ;

    options_description visible;
    visible.add(output).add(otc::standard_options());

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

std::ostream& operator<<(std::ostream& o, const RSplit& s) {
    writeSeparatedCollection(o, s.in, " ") <<" | ";
    if (s.out.size() < 100) {
        writeSeparatedCollection(o, s.out, " ");
    } else {
        auto it = s.out.begin();
        for(int i=0;i<100;i++) {
            o << *it++ <<" ";
        }
        o << "...";
    }
    return o;
}

/// Merge components c1 and c2 and return the component name that survived
int merge_components(int ic1, int ic2, vector<int>& component, vector<list<int>>& elements) {
    std::size_t c1 = static_cast<std::size_t>(ic1);
    std::size_t c2 = static_cast<std::size_t>(ic2);
    if (elements[c2].size() > elements[c1].size()) {
        std::swap(c1, c2);
    }
    for(int i: elements[c2]) {
        component[static_cast<std::size_t>(i)] = static_cast<int>(c1);
    }
    elements[c1].splice(elements[c1].end(), elements[c2]);
    return static_cast<int>(c1);
}

bool empty_intersection(const set<int>& xs, const vector<int>& ys) {
    for (int y: ys){
        if (xs.count(y)) {
            return false;
        }
    }
    return true;
}

static vector<int> indices;

/// Construct a tree with all the splits mentioned, and return a null pointer if this is not possible
unique_ptr<Tree_t> BUILD(const vector<int>& tips, const vector<const RSplit*>& splits) {
#pragma clang diagnostic ignored  "-Wsign-conversion"
#pragma clang diagnostic ignored  "-Wsign-compare"
#pragma clang diagnostic ignored  "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored  "-Wsign-compare"
    std::unique_ptr<Tree_t> tree(new Tree_t());
    tree->createRoot();
    // 1. First handle trees of size 1 and 2
    if (tips.size() == 1) {
        tree->getRoot()->setOttId(*tips.begin());
        return tree;
    } else if (tips.size() == 2) {
        auto Node1a = tree->createChild(tree->getRoot());
        auto Node1b = tree->createChild(tree->getRoot());
        auto it = tips.begin();
        Node1a->setOttId(*it++);
        Node1b->setOttId(*it++);
        return tree;
    }
    // 2. Initialize the mapping from elements to components
    vector<int> component;       // element index  -> component
    vector<list<int> > elements;  // component -> element indices
    for (int i=0;i<tips.size();i++) {
        indices[tips[i]] = i;
        component.push_back(i);
        elements.push_back({i});
    }
    // 3. For each split, all the leaves in the include group must be in the same component
    for(const auto& split: splits) {
        int c1 = -1;
        for(int i: split->in) {
            int j = indices[i];
            int c2 = component[j];
            if (c1 != -1 and c1 != c2) {
                merge_components(c1,c2,component,elements);
            }
            c1 = component[j];
        }
    }
    // 4. If we can't subdivide the leaves in any way, then the splits are not consistent, so return failure
    if (elements[component[0]].size() == tips.size()) {
        return {};
    }
    // 5. Make a vector of labels for the partition components
    vector<int> component_labels;                           // index -> component label
    vector<int> component_label_to_index(tips.size(),-1);   // component label -> index
    for (int c=0;c<tips.size();c++) {
        if (c == component[c]) {
            int index = component_labels.size();
            component_labels.push_back(c);
            component_label_to_index[c] = index;
        }
    }
    // 6. Create the vector of tips in each connected component 
    vector<vector<int>> subtips(component_labels.size());
    for(int i=0;i<component_labels.size();i++) {
        vector<int>& s = subtips[i];
        int c = component_labels[i];
        for (int j: elements[c]) {
            s.push_back(tips[j]);
        }
    }
    // 7. Determine the splits that are not satisfied yet and go into each component
    vector<vector<const RSplit*>> subsplits(component_labels.size());
    for(const auto& split: splits) {
        int first = indices[*split->in.begin()];
        assert(first >= 0);
        int c = component[first];
        // if none of the exclude group are in the component, then the split is satisfied by the top-level partition.
        bool satisfied = true;
        for(int x: split->out){
            if (indices[x] != -1 and component[indices[x]] == c) {
                satisfied = false;
                break;
            }
        }
        if (not satisfied) {
            int i = component_label_to_index[c];
            subsplits[i].push_back(split);
        }
    }
    // 8. Clear our map from id -> index, for use by subproblems.
    for(int id: tips) {
        indices[id] = -1;
    }
    // 9. Recursively solve the sub-problems of the partition components
    for(int i=0;i<subtips.size();i++) {
        auto subtree = BUILD(subtips[i], subsplits[i]);
        if (not subtree) {
            return {};
        }
        addSubtree(tree->getRoot(), *subtree);
    }
    return tree;
}

unique_ptr<Tree_t> BUILD(const vector<int>& tips, const vector<RSplit>& splits) {
    vector<const RSplit*> split_ptrs;
    for(const auto& split: splits) {
        split_ptrs.push_back(&split);
    }
    return BUILD(tips, split_ptrs);
}

/// Copy node names from taxonomy to tree based on ott ids, and copy the root name also
void add_names(unique_ptr<Tree_t>& tree, const unique_ptr<Tree_t>& taxonomy) {
    clearAndfillDesIdSets(*tree);
    for(auto n1: iter_post(*tree)) {
        for(auto n2: iter_post(*taxonomy)) {
            if (n1->getData().desIds == n2->getData().desIds) {
                n1->setName( n2->getName());
            }
        }
    }
}

set<int> remap_ids(const set<long>& s1, const map<long,int>& id_map) {
    set<int> s2;
    for(auto x: s1) {
        auto it = id_map.find(x);
        assert(it != id_map.end());
        s2.insert(it->second);
    }
    return s2;
}


/// Get the list of splits, and add them one at a time if they are consistent with previous splits
unique_ptr<Tree_t> combine(const vector<unique_ptr<Tree_t>>& trees, const set<long>& incertae_sedis, bool verbose) {
    // 0. Standardize names to 0..n-1 for this subproblem
    const auto& taxonomy = trees.back();
    auto all_leaves = taxonomy->getRoot()->getData().desIds;
    // index -> id
    vector<long> ids;
    // id -> index
    map<long,int> id_map;
    for(long id: all_leaves) {
        int i = ids.size();
        id_map[id] = i;
        ids.push_back(id);
        assert(id_map[ids[i]] == i);
        assert(ids[id_map[id]] == id);
    }
    auto remap = [&id_map](const set<long>& argIds) {return remap_ids(argIds, id_map);};
    vector<int> all_leaves_indices;
    for(int i=0;i<all_leaves.size();i++) {
        all_leaves_indices.push_back(i);
    }
    indices.resize(all_leaves.size());
    for(auto& i: indices) {
        i=-1;
    }
    // 1. Find splits in order of input trees
    vector<RSplit> consistent;
    for(const auto& tree: trees) {
        auto root = tree->getRoot();
        const auto leafTaxa = root->getData().desIds;
#pragma clang diagnostic ignored  "-Wunreachable-code-loop-increment"
        for(const auto& leaf: set_difference_as_set(leafTaxa, all_leaves)) {
            throw OTCError()<<"OTT Id "<<leaf<<" not in taxonomy!";
        }
        const auto leafTaxaIndices = remap(leafTaxa);
        for(auto nd: iter_post_const(*tree)) {
            if (nd->getData().desIds.size()>1 and leafTaxaIndices.size() > nd->getData().desIds.size()) {
                const auto descendants = remap(nd->getData().desIds);
                RSplit split{descendants, leafTaxaIndices};
                consistent.push_back(split);
                auto result = BUILD(all_leaves_indices, consistent);
                if (not result) {
                    consistent.pop_back();
                    if (verbose and nd->hasOttId()) {
                        LOG(INFO) << "Reject: ott" << nd->getOttId() << "\n";
                    }
                } else if (verbose and nd->hasOttId()) {
                    LOG(INFO) << "Keep: ott" << nd->getOttId() << "\n";
                }
            }
        }
    }
    // 2. Construct final tree and add names
    auto tree = BUILD(all_leaves_indices, consistent);
    for(auto nd: iter_pre(*tree)) {
        if (nd->isTip()) {
            int index = nd->getOttId();
            nd->setOttId(ids[index]);
        }
    }
    add_names(tree, taxonomy);
    return tree;
}

/// Create an unresolved taxonomy out of all the input trees.
unique_ptr<Tree_t> make_unresolved_tree(const vector<unique_ptr<Tree_t>>& trees, bool use_ids) {
    std::unique_ptr<Tree_t> retTree(new Tree_t());
    retTree->createRoot();
    if (use_ids) {
        map<long,string> names;
        for(const auto& tree: trees) {
            for(auto nd: iter_pre_const(*tree)) {
                if (nd->isTip()) {
                    long id = nd->getOttId();
                    auto it = names.find(id);
                    if (it == names.end()) {
                        names[id] = nd->getName();
                    }
                }
            }
        }
        for(const auto& n: names) {
            auto node = retTree->createChild(retTree->getRoot());
            node->setOttId(n.first);
            node->setName(n.second);
        }
        clearAndfillDesIdSets(*retTree);
    } else {
        set<string> names;
        for(const auto& tree: trees) {
            for(auto nd: iter_pre_const(*tree)) {
                if (nd->isTip()) {
                    names.insert(nd->getName());
                }
            }
        }
        for(const auto& n: names) {
            auto node = retTree->createChild(retTree->getRoot());
            node->setName(n);
        }
    }
    return retTree;
}

int main(int argc, char *argv[])
{
    try
    {
	// 1. Parse command line arguments
	variables_map args = parse_cmd_line(argc,argv);
  
	ParsingRules rules;
	rules.setOttIds = (bool)args.count("require-ott-ids");
	rules.pruneUnrecognizedInputTips = (bool)args.count("prune-unrecognized");

	bool synthesize_taxonomy = (bool)args.count("synthesize-taxonomy");
	bool cladeTips = not (bool)args.count("no-higher-tips");
	bool verbose = (bool)args.count("verbose");
	bool writeStandardized = (bool)args.count("standardize");
	if (writeStandardized) {
	    rules.setOttIds = false;
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

	//2.5 Load Incertae Sedis info
	std::set<long> incertae_sedis;
	if (args.count("incertae-sedis"))
	{
	    auto filename = args["incertae-sedis"].as<string>();
	    std::ifstream file(filename);
	    while (file)
	    {
		long i;
		file >> i;
		incertae_sedis.insert(i);
	    }
	}

	// 3. Make a fake taxonomy if asked
	if (synthesize_taxonomy) {
	    trees.push_back(make_unresolved_tree(trees, rules.setOttIds));
	    LOG(DEBUG)<<"taxonomy = "<<newick(*trees.back())<<"\n";
	}

	// 4. Add fake Ott Ids to tips and compute desIds (if asked)
	if (not rules.setOttIds) {
	    auto name_to_id = createIdsFromNames(*trees.back());
	    for(auto& tree: trees) {
		setIdsFromNamesAndRefresh(*tree, name_to_id);
	    }
	}

	// 5. Write out subproblem with newly minted ottids (if asked)
	if (writeStandardized) {
	    for(const auto& tree: trees) {
		relabelNodesWithOttId(*tree);
		std::cout<<newick(*tree)<<"\n";
	    }
	    return 0;
	}

	// 6. Check if trees are mapping to non-terminal taxa, and either fix the situation or die.
	for (int i = 0; i <trees.size() - 1; i++) {
	    if (cladeTips) {
		expandOTTInternalsWhichAreLeaves(*trees[i], *trees.back());
	    } else {
		requireTipsToBeMappedToTerminalTaxa(*trees[i], *trees.back());
	    }
	}

	// 7. Perform the synthesis
	auto tree = combine(trees, incertae_sedis, verbose);

	// 8. Set the root name (if asked)
	// FIXME: This could be avoided if the taxonomy tree in the subproblem always had a name for the root node.
	if (setRootName) {
	    tree->getRoot()->setName(args["root-name"].as<string>());
	}

	// 9. Write out the summary tree.
	writeTreeAsNewick(std::cout, *tree);
	std::cout<<"\n";

	return 0;
    }
    catch (std::exception& e)
    {
	std::cerr<<"otc-solve-subproblem: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}
