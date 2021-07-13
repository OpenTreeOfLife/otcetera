#include <algorithm>
#include <set>
#include <list>
#include <iterator>
#include <boost/smart_ptr/intrusive_ptr.hpp>

#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
#include "otc/tree_iter.h"
#include "otc/induced_tree.h"
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>
#include <optional>
#include <robin_hood.h>

#include "otc/conflict.h"

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

template <typename X>
using Set = robin_hood::unordered_set<X>;

using namespace otc;

typedef TreeMappedWithSplits Tree_t;
typedef Tree_t::node_type node_t;

static vector<int> indices;

int depth(const Tree_t::node_type* nd)
{
    return nd->get_data().depth;
}

/// Create a SORTED vector from a set
template <typename T>
vector<T> set_to_vector(const set<T>& s) {
    vector<T> v;
    v.reserve(s.size());
    std::copy(s.begin(), s.end(), std::back_inserter(v));
    return v;
}

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

struct RSplitObj
{
    static std::size_t num;

    mutable int _refs = 0;

    friend inline void intrusive_ptr_release(RSplitObj* pThis)
    {
        if (--pThis->_refs == 0 ) {
            delete pThis;
        }
    }

    friend inline void intrusive_ptr_add_ref(RSplitObj* pThis)
    {
        pThis->_refs++;
    }

    friend inline void intrusive_ptr_release(const RSplitObj* pThis)
    {
        if(--const_cast<RSplitObj*>(pThis)->_refs == 0 ) {
            delete const_cast<RSplitObj*>(pThis);
        }
    }

    friend inline void intrusive_ptr_add_ref(const RSplitObj* pThis)
    {
        const_cast<RSplitObj*>(pThis)->_refs++;
    }

    vector<int> in;
    vector<int> out;
    optional<std::size_t> id;

    RSplitObj()
    {
        id = num++;
    }

    RSplitObj(const set<int>& i, const set<int>& a)
    {
        id = num++;
        in  = set_to_vector(i);
        set_difference(begin(a), end(a), begin(in), end(in), std::inserter(out, out.end()));
    }
};

std::size_t RSplitObj::num = 0;

using RSplit = boost::intrusive_ptr<RSplitObj>;
using ConstRSplit = boost::intrusive_ptr<const RSplitObj>;

std::ostream& operator<<(std::ostream& o, const RSplitObj& s)
{
    o<<"["<<s.in.size() + s.out.size()<<" tips] ";
    for(auto& is: s.in)
        o<<"ott"<<is<<" ";
    o<<"| ";
    for(auto& os: s.out)
        o<<"ott"<<os<<" ";
    return o;
}

RSplit split_from_include_exclude(const set<int>& i, const set<int>& e)
{
    RSplit s(new RSplitObj);
    s->in = set_to_vector(i);
    s->out = set_to_vector(e);
    return s;
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
        ("incertae-sedis,I",value<string>(),"File containing Incertae sedis ids")
        ("root-name,n",value<string>(), "Rename the root to this name")
        ("no-higher-tips,l", "Tips may be internal nodes on the taxonomy.")
        ("prune-unrecognized,p","Prune unrecognized tips");
    
    options_description strategies("Solver strategies");
    strategies.add_options()
        ("batching",value<bool>()->default_value(true), "Make unresolved taxonomy from input tips.")
        ("oracle", value<bool>()->default_value(true), "Predict conflicting splits before BUILD.")
        ("incremental", value<bool>()->default_value(true),"Reuse work from previous BUILD.")
        ;

    options_description other("Other options");
    other.add_options()
        ("synthesize-taxonomy,T","Make unresolved taxonomy from input tips.")
        ("allow-no-ids,a", "Allow problems w/o OTT ids")
        ("standardize,S", "Write out a standardized subproblem and exit.")
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

std::ostream& operator<<(std::ostream& o, const ConstRSplit& s) {
    write_separated_collection(o, s->in, " ") <<" | ";
    if (s->out.size() < 100) {
        write_separated_collection(o, s->out, " ");
    } else {
        auto it = s->out.begin();
        for(int i=0;i<100;i++) {
            o << *it++ <<" ";
        }
        o << "...";
    }
    return o;
}

std::ostream& operator<<(std::ostream& o, const RSplit& s) {
    return o<<(ConstRSplit(s));
}

struct Solution;

struct component_t
{
    list<int> elements;

    shared_ptr<Solution> solution;
    vector<shared_ptr<Solution>> old_solutions;

    vector<int> get_taxa(const std::vector<int>& other_taxa) const
    {
        vector<int> taxa;
        taxa.reserve(elements.size());
        for(auto index: elements)
            taxa.push_back(other_taxa[index]);
        return taxa;
    }
};

template <typename T>
bool sort_cmp(const vector<T>& v1,  const vector<T>& v2)
{
    if (v1.size() != v2.size()) return false;

    auto w1 = v1; std::sort(w1.begin(), w1.end());
    auto w2 = v2; std::sort(w2.begin(), w2.end());

    return (w1 == w2);
}

typedef component_t* component_ref;

// A "partial" solution only has implied_splits + non_implied_splits
// A "full" solution also has components.
// Eventually I would like to replace non_implied_splits with the components!

struct Solution
{
    vector<int> taxa;

    vector<ConstRSplit> implied_splits;      // alpha
    vector<ConstRSplit> non_implied_splits;  // beta

    vector< component_ref > component_for_index;
    vector< unique_ptr<component_t> > components;

    vector<ConstRSplit> non_implied_splits_from_components() const;
    vector<ConstRSplit> splits_from_components() const;

    unique_ptr<Tree_t> get_tree() const;

    Solution(const vector<int>& t, const vector<ConstRSplit>& s)
        :taxa(t), non_implied_splits(s), component_for_index(taxa.size())
    {}
    Solution(const vector<int>& t, vector<ConstRSplit>&& s)
        :taxa(t), non_implied_splits(std::move(s)), component_for_index(taxa.size())
    {}
    Solution(const vector<int>& t)
        :Solution(t,{})
    {}
    Solution(const component_t& c, const std::vector<int> other_taxa)
        :Solution(c.get_taxa(other_taxa))
    {}
};

vector<ConstRSplit> Solution::splits_from_components() const
{
    vector<ConstRSplit> splits = non_implied_splits_from_components();
    for(auto& split: implied_splits)
        splits.push_back(split);
    return splits;
}

vector<ConstRSplit> Solution::non_implied_splits_from_components() const
{
    vector<ConstRSplit> splits;
    for(auto& component: components)
    {
        for(auto split: component->solution->splits_from_components())
            splits.push_back(split);
    }
    return splits;
}

/// Merge components c1 and c2 and return the component name that survived
void merge_component_with_trivial(component_ref c1, int index2, vector<component_ref>& component)
{
    component[index2] = c1;
    c1->elements.push_back(index2);


    c1->solution = {};
}

bool exclude_group_intersects_component(const ConstRSplit& split, const component_t* component, const vector<component_ref>& component_for_index)
{
    for(int taxon: split->out)
    {
        int index = indices[taxon];

        if (index == -1) continue;

        if (component_for_index[index] == component) return true;
    }
    return false;
}

template <typename T>
void append(vector<T>& v1, const vector<T>& v2)
{
    v1.insert(v1.end(), v2.begin(), v2.end());
}

/// Merge components c1 and c2 and return the component name that survived
component_ref merge_components(component_ref c1, component_ref c2, vector<component_ref>& component)
{
    if (c2->elements.size() > c1->elements.size())
        std::swap(c1, c2);

    for(int i: c2->elements)
        component[i] = c1;

    c1->elements.splice(c1->elements.end(), c2->elements);

    // One of these components could be new -- that is, composed only of previously-trivial components.
    append(c1->old_solutions, c2->old_solutions);

    c1->solution = {};

    return c1;
}

bool empty_intersection(const set<int>& xs, const vector<int>& ys) {
    for (int y: ys){
        if (xs.count(y)) {
            return false;
        }
    }
    return true;
}

template <typename T>
T remove_unordered(std::vector<T>& v, int i)
{
    assert(0 <= i and i < v.size());

    auto t = v[i];
    if (i < int(v.size())-1)
        std::swap(v[i], v.back());
    v.pop_back();

    return t;
}

unique_ptr<Tree_t> Solution::get_tree() const
{
    assert(taxa.size() > 1);

    // 1. Make a tree with just a root node
    std::unique_ptr<Tree_t> tree(new Tree_t());
    tree->create_root();

    // 2. Add children for non-trivial components
    for(auto& component: components)
        add_subtree(tree->get_root(), *component->solution->get_tree());

    // 3. Add children for trivial components
    for(int index=0;index<taxa.size();index++)
    {
        if (not component_for_index[index])
        {
            auto taxon = taxa[index];
            auto node = tree->create_child(tree->get_root());
            node->set_ott_id(taxon);
        }
    }

    return tree;
}

// 0. Update doc?

// 1. Handle adding trivial taxa to a component.

// 2. When merging components, we need to consider moving splits from implied -> non_implied.
// Should we do this from scratch?

// 3. Each new component should have a list of pointers to previous .. components? solutions?
// We could get their splits/taxa all at once, instead of copying them multiple times.


/// Construct a tree with all the splits mentioned, and return a null pointer if this is not possible
bool BUILD(Solution& solution)
{
#pragma clang diagnostic ignored  "-Wsign-conversion"
#pragma clang diagnostic ignored  "-Wsign-compare"
#pragma clang diagnostic ignored  "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored  "-Wsign-compare"

    // This describes the problem
    auto& taxa = solution.taxa;
    auto& new_splits = solution.non_implied_splits;

    auto& component_for_index = solution.component_for_index;
    auto& components = solution.components;

    // 1. If there are no splits, then we are consistent.
    if (new_splits.empty())
        return true;

    // 2. Initialize the mapping from taxa to indices.
    for(int k=0;k<indices.size();k++)
        assert(indices[k] == -1);
    for (int i=0;i<taxa.size();i++)
        indices[taxa[i]] = i;

    // 3. For each split, all the leaves in the include group must be in the same component
    for(const auto& split: new_splits)
    {
        component_ref split_comp = nullptr;
        for(int taxon: split->in)
        {
            int index = indices[taxon];
            assert(index != -1);
            auto taxon_comp = component_for_index[index];
            if (not split_comp)
            {
                if (not taxon_comp)
                {
                    components.push_back(std::make_unique<component_t>());
                    taxon_comp = components.back().get();
                    merge_component_with_trivial(taxon_comp, index, component_for_index);
                }
                split_comp = taxon_comp;
            }
            else if (not taxon_comp)
                merge_component_with_trivial(split_comp, index, component_for_index);
            else if (split_comp != taxon_comp)
                split_comp = merge_components(split_comp,taxon_comp,component_for_index);
        }
    }

    // 4. If we can't subdivide the leaves in any way, then the splits are not consistent, so return failure
    if (component_for_index[0] and component_for_index[0]->elements.size() == taxa.size())
    {
        for(int id: taxa)
            indices[id] = -1;
        return false;
    }

    // 5. Pack the components
    vector<unique_ptr<component_t>> packed_components;
    for(auto& component: components)
        if (not component->elements.empty())
            packed_components.push_back( std::move(component) );
    std::swap(components, packed_components);


    // OK, so can we change how component merges work?
    // When do components first get (i) old_implied_splits, (ii) old_non_implied_splits, and (iii) new_splits?

    // I would like to change how component merges work.
    // - If we merge two components that already have solutions, then each solution/old component should keep its implied_splits + non_implied_splits.
    // - When handling each component, we want to look at (a) new splits to the component and (b) splits for merged components.
    // - So, maybe the implied_splits + non_implied_splits should be part of the solution object...

    // 6. Check implied splits to see if they are STILL implied.
    for(auto& component: components)
    {
        // We don't need to re-check implied_splits if the taxon set hasn't changed.
        if (component->solution)
        {
            assert((component->old_solutions.size() == 1) and (component->elements.size() == component->old_solutions[0]->taxa.size()));
            continue;
        }

        assert(not component->solution);
        component->solution = std::make_shared<Solution>(*component, taxa);

        // FIXME: instead of collecting implied splits and non-implied splits,
        //        only scan the implied splits to see if the component has been punctured.
        //        If the component is punctured, then we divide the A[i] between A' and B', and
        //          add the child solutions to our list.
        //        We should also be able to GENERATE the non_implied splits B[i] using old_solution->non_implied_splits_from_components()
        for(auto& old_solution: component->old_solutions)
        {
            append(component->solution->implied_splits, old_solution->implied_splits);
            auto non_implied_splits = old_solution->non_implied_splits_from_components();
            append(component->solution->non_implied_splits, non_implied_splits);
        }

        for(int i=0;i<component->solution->implied_splits.size();)
        {
            auto& split = component->solution->implied_splits[i];
#ifndef NDEBUG
            int first = indices[*split->in.begin()];
            assert(first >= 0);
            assert(component_for_index[first] == component.get());
#endif
            bool implied = not exclude_group_intersects_component(split, component.get(), component_for_index);
            if (not implied)
            {
                auto split = remove_unordered(component->solution->implied_splits,i);
                component->solution->non_implied_splits.push_back(split);
            }
            else
                i++;
        }
    }

    // NOTE: If all new splits are implied and no OLD splits are implied, then perhaps
    //       all the old components will be sub-problems of the merged component?

    // 7. Determine the splits that are not satisfied yet and go into each component
    for(auto& split: new_splits)
    {
        int first = indices[*split->in.begin()];
        assert(first >= 0);
        auto component = component_for_index[first];

        bool implied = not exclude_group_intersects_component(split, component, component_for_index);
        if (implied)
            component->solution->implied_splits.push_back(split);
        else
            component->solution->non_implied_splits.push_back(split);
    }

    // 8. We've now set up the sub-problems, so we can clear non_implied_splits/new_splits.
    solution.non_implied_splits.clear();

    // 9. Clear our map from id -> index, for use by subproblems.
    for(int id: taxa) {
        indices[id] = -1;
    }

    // 10. Recursively solve the sub-problems of the partition components
    for(auto& component: components)
    {
        assert(component->elements.size() >= 2);

        if (not BUILD(*component->solution))
            return false;

        component->old_solutions = { component->solution };
    }

    return true;
}

bool BUILD(Solution& solution, const vector<ConstRSplit>& new_splits)
{
    for(auto split: new_splits)
        solution.non_implied_splits.push_back(split);

    return BUILD(solution);
}

template <typename T>
bool is_subset(const std::set<T>& set_two, const std::set<T>& set_one) {
    return std::includes(set_one.begin(), set_one.end(), set_two.begin(), set_two.end());
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

set<int> remap_ids(const set<OttId>& s1, const map<OttId,int>& id_map) {
    set<int> s2;
    for(auto x: s1) {
        auto it = id_map.find(x);
        assert(it != id_map.end());
        s2.insert(it->second);
    }
    return s2;
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


set<Tree_t::node_type*> find_conflicting_nodes(unique_ptr<Tree_t>& ok_tree, unique_ptr<Tree_t>& tree_to_clean)
{
    set<node_t*> conflicting_nodes;

    typedef Tree_t::node_type node_t;
    typedef ConflictTree::node_type cnode_t;
    std::function<node_t*(node_t*,node_t*)> mrca_of_pair = [](node_t* n1, node_t* n2) {return mrca_from_depth(n1,n2);};

    auto tree_to_clean_ottid_to_node = get_ottid_to_node_map(*tree_to_clean);

    auto ok_tree_ottid_to_node = get_ottid_to_node_map(*ok_tree);

    compute_depth(*ok_tree);

    compute_depth(*tree_to_clean);

    auto [induced_tree_to_clean,to_induced] = get_induced_tree_and_node_map<ConflictTree>(*tree_to_clean,
                                                                                          tree_to_clean_ottid_to_node,
                                                                                          mrca_of_pair,
                                                                                          *ok_tree,
                                                                                          ok_tree_ottid_to_node);

    auto induced_ok_tree = get_induced_tree<ConflictTree>(*ok_tree,
                                                          ok_tree_ottid_to_node,
                                                          mrca_of_pair,
                                                          *tree_to_clean,
                                                          tree_to_clean_ottid_to_node);

    if (induced_tree_to_clean->get_root() and induced_ok_tree->get_root())
    {
        std::unordered_map<const cnode_t*, node_t*> from_induced;
        for(auto& [non_induced,induced]: to_induced)
        {
            from_induced[induced] = non_induced;
        }

        auto log_supported_by    = [&](const cnode_t* /* node2 */, const cnode_t* /* node1 */) {};
        auto log_partial_path_of = [&](const cnode_t* /* node2 */, const cnode_t* /* node1 */) {};
        auto log_resolved_by     = [&](const cnode_t* /* node2 */, const cnode_t* /* node1 */) {};
        auto log_terminal        = [&](const cnode_t* /* node2 */, const cnode_t* /* node1 */) {};
        auto log_conflicts_with  = [&](const cnode_t* /* node2 */, const cnode_t* node1 )
        {
            assert(node1->has_children());
            // for each of node1 and all its monotypic ancestors.
            do
            {
                // mark the corresponding node in the original (not induced) tree as conflicting.
                auto node2 = from_induced.at(node1);
                assert(node2->has_children());
                conflicting_nodes.insert(node2);

                node1 = node1->get_parent();
            }
            while(node1->is_outdegree_one_node());
        };

        perform_conflict_analysis(*induced_tree_to_clean, *induced_ok_tree, log_supported_by, log_partial_path_of, log_conflicts_with, log_resolved_by, log_terminal);
    }

    return conflicting_nodes;
}

template<typename N>
inline void collapse_node_(N* nd)
{
    assert(nd);
    assert(nd->has_children());

    while(nd->has_children())
    {
        auto child = nd->get_first_child();
        child->detach_this_node();
        nd->add_sib_on_left(child);
    }
    nd->detach_this_node();
}


void remove_conflicting_splits_from_tree(unique_ptr<Tree_t>& ok_tree, unique_ptr<Tree_t>& tree_to_clean)
{
    for(auto& node: find_conflicting_nodes(ok_tree, tree_to_clean))
        collapse_node_(node);
}

void remove_conflicting_splits_from_tree(vector<unique_ptr<Tree_t>>& trees, int k)
{
    for(int i=0;i<k;i++)
        remove_conflicting_splits_from_tree(trees[i],trees[k]);
}

bool conflicting(const vector<int>& all_leaves_indices, const vector<ConstRSplit>& splits)
{
    auto solution = std::make_shared<Solution>(all_leaves_indices);
    auto result = BUILD(*solution, splits);
    return not result;
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

/// Get the list of splits, and add them one at a time if they are consistent with previous splits
unique_ptr<Tree_t> combine(vector<unique_ptr<Tree_t>>& trees, const set<OttId>& incertae_sedis, variables_map& args)
{
    bool verbose = (bool)args.count("verbose");
    bool batching = args["batching"].as<bool>();
    bool oracle = args["oracle"].as<bool>();
    bool incremental = args["incremental"].as<bool>();

    // 1. Standardize names to 0..n-1 for this subproblem
    const auto& taxonomy = trees.back();
    auto all_leaves = taxonomy->get_root()->get_data().des_ids;

    // ids:    index -> id
    // id_map: id    -> index
    auto [ids, id_map] = make_index_map(all_leaves);

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
            if (incremental and solution)
            {
                vector<ConstRSplit> new_splits;
                for(int i=0;i<n;i++)
                    new_splits.push_back(splits[start+i].second);

                result = BUILD(*solution, new_splits);

                if (result)
                {
                    for(auto& new_split: new_splits)
                        consistent.push_back(new_split);
                }
            }
            else
            {
                for(int i=0;i<n;i++)
                    consistent.push_back(splits[start+i].second);

                solution = std::make_shared<Solution>(all_leaves_indices, consistent);

                result = BUILD(*solution);

                if (not result)
                {
                    for(int i=0;i<n;i++)
                        consistent.pop_back();
                }
            }

            total_build_calls ++;

            if (not incremental or not result) solution = {};

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
            for(int i=0;i<splits2.size();i++)
                add_splits_if_consistent_batch(splits2,i,1);
        }

        LOG(INFO)<<"i = "<<i<<"  Total build calls = "<<total_build_calls;
    }

    vector<const_node_type<Tree_t>*> compatible_taxa;
    for(auto node: iter_pre(*trees.back()))
    {
        if (node->get_parent() and not node->is_tip())
            compatible_taxa.push_back(node);
    }

    // 2. Construct final tree and add names

    //FIXME - discard previous solution;
    solution = std::make_shared<Solution>(all_leaves_indices, consistent);
    auto result = BUILD(*solution);
    assert(result);
    auto tree = solution->get_tree();
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
    return tree;
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
        write_tree_as_newick(std::cout, *tree);
        std::cout << "\n";
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
