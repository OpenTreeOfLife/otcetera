#include "combine.h"

#include "build.h"
#include "oracle.h"

#include <chrono>
#include <map>
#include <utility>
#include <iomanip>

using namespace otc;

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

using boost::program_options::variables_map;

Tree_t::node_type* find_mrca_of_desids(const set<OttId>& ids, const std::unordered_map<OttId, Tree_t::node_type*>& summaryOttIdToNode) {
    int first = *ids.begin();
    auto node = summaryOttIdToNode.at(first);
    while( not is_subset(ids, node->get_data().des_ids) )
        node = node->get_parent();
    return node;
}

bool is_ancestor_of(const Tree_t::node_type* n1, const Tree_t::node_type* n2) {
    // make sure the depth fields are initialized
    assert(n1 == n2 or depth(n1) != 0 or depth(n2) != 0);
    if (depth(n2) > depth(n1)) {
        while (depth(n2) != depth(n1)) {
            n2 = n2->get_parent();
        }
        return (n2 == n1);
    } else {
        return false;
    }
}


const Tree_t::node_type* find_unique_maximum(const vector<const Tree_t::node_type*>& nodes)
{
    for(int i=0;i<nodes.size();i++)
    {
        bool is_ancestor = true;
        for(int j=0;j<nodes.size() and is_ancestor;j++) {
            if (j==i) {
                continue;
            }
            if (not is_ancestor_of(nodes[i],nodes[j])) {
                is_ancestor = false;
            }
        }
        if (is_ancestor) {
            return nodes[i];
        }
    }
    return nullptr;
}

Tree_t::node_type* add_monotypic_parent(Tree_t& tree, Tree_t::node_type* nd) {
    if (nd->get_parent()) {
        auto p = nd->get_parent();
        auto monotypic = tree.create_child(p);
        nd->detach_this_node();
        monotypic->add_child(nd);
        return monotypic;
    } else {
        auto monotypic = tree.create_root();
        monotypic->add_child(nd);
        return monotypic;
    }
}

const Tree_t::node_type* select_canonical_ottid(const vector<const Tree_t::node_type*>& names) {
    // We should only have to make this choice if there are at least 2 names to choose from.
    assert(names.size() >= 2);
    // Do something more intelligent here - perhaps prefer non-incertae-sedis taxa, and then choose lowest ottid.
    return names.front();
}

void register_ottid_equivalences(const Tree_t::node_type* canonical, const vector<const Tree_t::node_type*>& names) {
    // First pass - actually we should write on a JSON file.
    std::cerr << canonical->get_name() << " (canonical): equivalent to ";
    for(auto name: names)
        std::cerr << name->get_name() << " ";
    std::cerr << "\n";
}

void add_root_and_tip_names(Tree_t& summary, Tree_t& taxonomy) {
    // name root
    summary.get_root()->set_name(taxonomy.get_root()->get_name());
    if (taxonomy.get_root()->has_ott_id()) {
        summary.get_root()->set_ott_id(taxonomy.get_root()->get_ott_id());
    }
    // name tips
    auto summaryOttIdToNode = get_ottid_to_node_map(summary);
    for(auto nd: iter_leaf(taxonomy)) {
        auto id = nd->get_ott_id();
        auto nd2 = summaryOttIdToNode.at(id);
        nd2->set_name( nd->get_name());
    }
}

//  Given a list of taxon nodes that are known NOT to conflict with the summary tree,
//   * map each name to the MRCA of its include group.
//   * copy the (i) node name and (ii) ottid to the MRCA on the summary tree.
//
//  If multiple names get assigned to a single node, try to handle this by
//    making a monotypic parent assigning the rootmost name to that.
//  Otherwise choose a name arbitrarily.
void add_names(Tree_t& summary, const vector<const Tree_t::node_type*>& compatible_taxa) {
    auto summaryOttIdToNode = get_ottid_to_node_map(summary);

    // 1. Set the des_ids for the summary
    clear_and_fill_des_ids(summary);
    // 2. Place each taxon N at the MRCA of its include group.
    map<Tree_t::node_type*, vector<const Tree_t::node_type*>> name_groups;
    for(auto n2: compatible_taxa) {
        auto mrca = find_mrca_of_desids(n2->get_data().des_ids, summaryOttIdToNode);
        if (not name_groups.count(mrca)) {
            name_groups[mrca] = {};
        }
        name_groups[mrca].push_back(n2);

        // Any extra desids are here because an incertae sedis taxon was placed inside this node,
        // or inside a child.
    }
    // 3. Handle each summary 
    for(auto& name_group: name_groups) {
        auto summary_node = name_group.first;
        auto & names = name_group.second;
        // 3.1. As long as there is a unique root-most name, put that name in a monotypic parent.
        // This can occur when a node has two children, and one of them is an incertae sedis taxon that is moved more tip-ward.
        while (auto max = find_unique_maximum(names)) {
            if (names.size() == 1) {
                set_name_and_maybe_ott_id(*max, *summary_node);
            } else {
                auto p = add_monotypic_parent(summary, summary_node);
                set_name_and_maybe_ott_id(*max, *p);
                p->get_data().des_ids = p->get_first_child()->get_data().des_ids;
            }
            // Move the "removed" elements to the end and the erase them.  Weird.
            names.erase(std::remove(names.begin(), names.end(), max), names.end());
        }
        // 3.2. Select a canonical name from the remaining names.
        if (not names.empty()) {
            // Select a specific ottid as the canonical name for this summary node
            auto canonical = select_canonical_ottid(names);
            summary_node->set_name(canonical->get_name());

            // Write out the equivalence of the remaining ottids to the canonical ottid
            names.erase(std::remove(names.begin(), names.end(), canonical), names.end());
            register_ottid_equivalences(canonical, names);
        }
    }
}

template <typename Tree_t>
vector<typename Tree_t::node_type const*> get_siblings(typename Tree_t::node_type const* nd) {
    vector<typename Tree_t::node_type const*> sibs;
    for(auto sib = nd->get_first_sib(); sib; sib = sib->get_next_sib()) {
        if (sib != nd) {
            sibs.push_back(sib);
        }
    }
    return sibs;
}

template<typename Tree_T>
map<typename Tree_t::node_type const*, set<OttId>> construct_include_sets(const Tree_t& tree, const set<OttId>& incertae_sedis)
{
    map<typename Tree_t::node_type const*, set<OttId>> include;
    for(auto nd: iter_post(tree)) {
        // 1. Initialize set for this node.
        auto & inc = include[nd];
        // 2. Add OttId for tip nodes
        if (nd->is_tip()) {
            inc.insert(nd->get_ott_id());
        } else if (nd == tree.get_root()) {
            continue;
        }
        // 3. Add Ids of children only if they are NOT incertae sedis
        for(auto nd2: iter_child(*nd)) {
            if (not incertae_sedis.count(nd2->get_ott_id())) {
                auto& inc_child = nd2->get_data().des_ids;
                inc.insert(begin(inc_child),end(inc_child));
            }
        }
    }
    return include;
}

template<typename Tree_t>
map<typename Tree_t::node_type const*, set<OttId>> construct_exclude_sets(const Tree_t& tree, const set<OttId>& incertae_sedis) {
    map<typename Tree_t::node_type const*, set<OttId>> exclude;
    // 1. Set exclude set for root node to the empty set.
    exclude[tree.get_root()];
    for(auto nd: iter_pre(tree)) {
        // 2. Skip tips and the root node.
        if (nd->is_tip() || nd == tree.get_root()) {
            continue;
        }
        // 3. Start with the exclude set for the parent.  This should already exist.
        set<OttId> ex = exclude.at(nd->get_parent());
        // 4. The exclude set should ALSO include ALL (not just some) descendants of siblings.
        for(auto nd2: get_siblings<Tree_t>(nd)) {
            if (not incertae_sedis.count(nd2->get_ott_id())) {
                // 5. In this variant, we DO exclude descendants that are accessed through a node marked I.S.
                auto& ex_sib = nd2->get_data().des_ids;
                ex.insert(begin(ex_sib),end(ex_sib));
            }
        }
        exclude[nd] = ex;
    }
    return exclude;
}

vector<pair<node_type<Tree_t>*,RSplit>>
splits_for_tree(Tree_t& tree, const std::function< set<int>(const set<OttId>&) >& remap)
{
    vector<pair<node_type<Tree_t>*,RSplit>> splits;
    auto root = tree.get_root();
    const auto leafTaxa = root->get_data().des_ids;
    const auto leafTaxaIndices = remap(leafTaxa);
    for(auto nd: iter_pre(tree))
    {
        if (not nd->is_tip() and nd != root)
        {
            auto descendants = remap(nd->get_data().des_ids);
            splits.push_back({nd,RSplit(new RSplitObj{descendants, leafTaxaIndices})});
        }
    }
    return splits;
}

vector<pair<node_type<Tree_t>*,RSplit>>
splits_for_taxonomy_tree(Tree_t& tree, const std::function< set<int>(const set<OttId>&) >& remap, const set<OttId>& incertae_sedis)
{
    if (incertae_sedis.empty())
        return splits_for_tree(tree, remap);

    vector<pair<node_type<Tree_t>*,RSplit>> splits;
    auto root = tree.get_root();

    auto exclude = construct_exclude_sets(tree, incertae_sedis);

    for(auto nd: iter_pre(tree))
    {
        if (not nd->is_tip() and nd != root) {
            // construct split
            const auto descendants = remap(nd->get_data().des_ids);
            const auto nondescendants = remap(exclude[nd]);
            splits.push_back({nd, split_from_include_exclude(descendants, nondescendants)});
        }
    }

    return splits;
}

set<int> remap_ids(const set<OttId>& s1, const map<OttId,int>& id_map) {
    set<int> s2;
    for(auto x: s1) {
        auto it = id_map.find(x);
        assert(it != id_map.end());
        s2.insert(it->second);
    }
    return s2;
}

template <typename T>
pair<vector<T>,map<T,int>> make_index_map(const set<T>& s)
{
    pair<vector<T>,map<T,int>> x;
    // ids:    index -> id
    // id_map: id    -> index
    auto& [ids,id_map] = x;
    for(auto& id: s)
    {
        int i = ids.size();
        id_map[id] = i;
        ids.push_back(id);
        assert(id_map[ids[i]] == i);
        assert(ids[id_map[id]] == id);
    }
    return x;
}


/// Get the list of splits, and add them one at a time if they are consistent with previous splits
unique_ptr<Tree_t> combine(vector<unique_ptr<Tree_t>>& trees, const set<OttId>& incertae_sedis, variables_map& args)
{
    bool verbose = (bool)args.count("verbose");
    bool batching = args["batching"].as<bool>();
    bool oracle = args["oracle"].as<bool>();
    bool incremental = args["incremental"].as<bool>();
    bool g_do_timing = (bool)args.count("time");
    auto start_timing = std::chrono::high_resolution_clock::now();
        
    // 1. Standardize names to 0..n-1 for this subproblem
    const auto& taxonomy = trees.back();
    auto all_leaves = taxonomy->get_root()->get_data().des_ids;

    // ids:    index -> id
    // id_map: id    -> index
    auto [ids, id_map_] = make_index_map(all_leaves);
    auto& id_map = id_map_;

    std::function< set<int>(const set<OttId>&) > remap = [&id_map](const set<OttId>& argIds) {return remap_ids(argIds, id_map);};
    vector<int> all_leaves_indices;
    for(int i=0;i<all_leaves.size();i++) {
        all_leaves_indices.push_back(i);
    }
    indices.resize(all_leaves.size());
    for(auto& i: indices) {
        i=-1;
    }

    /// Incrementally add splits from @splits_to_try to @consistent if they are consistent with it.
    vector<ConstRSplit> consistent;

    shared_ptr<Solution> solution;

    // A non-null solution means that consistent splits are already part of the solution.

    int total_build_calls = 0;
    auto add_splits_if_consistent = [&](vector<pair<node_type<Tree_t>*,RSplit>>& splits, int start, int n)
        {
            bool result;
            if (incremental)
            {
                if (not solution)
                    solution = std::make_shared<Solution>(all_leaves_indices);
                vector<ConstRSplit> new_splits;
                for(int i=0;i<n;i++)
                    new_splits.push_back(splits[start+i].second);

                result = BUILDINC(solution, new_splits);
                LOG(TRACE)<<"consistent = "<< consistent.size()<<" -> "<<consistent.size()+n<<": "<<(result?"ok":"FAIL");
                if (result)
                {
                    for(auto& new_split: new_splits)
                        consistent.push_back(new_split);
                }
                else
                {
                    LOG(TRACE)<<"FAIL!";
                }
                assert(consistent.size() == solution->n_splits_from_components());
            }
            else
            {
                for(int i=0;i<n;i++)
                    consistent.push_back(splits[start+i].second);

                solution = std::make_shared<Solution>(all_leaves_indices);

                result = BUILDINC(solution, consistent);
                LOG(TRACE)<<"consistent = "<< consistent.size()-n<<" -> "<<consistent.size()<<": "<<(result?"ok":"FAIL");
                if (not result)
                {
                    for(int i=0;i<n;i++)
                        consistent.pop_back();
                }
            }

            total_build_calls ++;

            if (n==1 and not result) collapse_node_(splits[start].first);

            return result;
        };

    std::function<void(vector<pair<node_type<Tree_t>*,RSplit>>&,int,int)> add_splits_if_consistent_batch;
    add_splits_if_consistent_batch = [&](vector<pair<node_type<Tree_t>*,RSplit>>& splits, int start, int n)
        {
            assert(n >= 1);
            assert(start+n <= splits.size());
            auto result = add_splits_if_consistent(splits, start, n);
            if (not result and n > 1)
            {
                int n1 = n/2;
                int n2 = n - n1;
                add_splits_if_consistent_batch(splits, start   , n1);
                add_splits_if_consistent_batch(splits, start+n1, n2);
            }
        };

    // 1. Find splits in order of input trees
    vector<pair<node_type<Tree_t>*,RSplit>> splits;
    for(int i=0;i<trees.size();i++)
    {
        const auto& tree = trees[i];

        // 1. Remove splits from tree i that directly conflict with previous TREES.
        //    Unless this is the taxonomy tree and there are incertae sedis taxa.
        if (oracle and (i<int(trees.size())-1 or incertae_sedis.empty()))
            remove_conflicting_splits_from_tree(trees,i);

        // 2. Get remaining splits
        auto splits2 = (i<trees.size()-1)
            ?splits_for_tree(*tree, remap)
            :splits_for_taxonomy_tree(*tree, remap, incertae_sedis);

        if (splits2.empty()) continue;

        // 3. Add compatible splits to `splits` and remove incompatible nodes from `trees[i]`;
        if (batching)
            add_splits_if_consistent_batch(splits2, 0, splits2.size());
        else
        {
            for(int j=0;j<splits2.size();j++)
                add_splits_if_consistent_batch(splits2,j,1);
        }

        LOG(DEBUG)<<"i = "<<i<<"  Total build calls = "<<total_build_calls;
    }

    vector<const_node_type<Tree_t>*> compatible_taxa;
    for(auto node: iter_pre(*trees.back()))
    {
        if (node->get_parent() and not node->is_tip())
            compatible_taxa.push_back(node);
    }

    // 2. Construct final tree and add names

    //FIXME - discard previous solution;
    auto tree = BUILD(all_leaves_indices, consistent);
    assert(tree);

    for(auto nd: iter_pre(*tree))
    {
        if (nd->is_tip())
        {
            int index = nd->get_ott_id();
            nd->set_ott_id(ids[index]);
        }
    }

    // We've modified the local copy of the taxonomy, so recompute depths.
    // "compatible_taxa" has pointers into it, but should only have pointers to surviving nodes.
    compute_depth(*taxonomy);
    add_root_and_tip_names(*tree, *taxonomy);
    add_names(*tree, compatible_taxa);
    if (g_do_timing) {
        auto end_timing = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end_timing - start_timing;
        std::cerr << "timing " << std::setprecision(10) << diff.count() << " seconds.\n";
    }
    return tree;
}

