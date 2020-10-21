#include <algorithm>
#include <set>
#include <list>
#include <iterator>

#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
#include "otc/tree_iter.h"
#include "otc/induced_tree.h"
#include <fstream>
#include <queue>
#include <sstream>
#include <boost/filesystem.hpp>
#include "robin_hood.h"
#include <random>
#include <limits>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <optional>

/* TODO:
 *
 * Correctness:
 * C1. We need to walk the tree from the root toward the tips,
 *       and after we clean out a subtree, include it in the list when
 *       we check sibling subtrees.
 *
 * Speed:
 * S1. We need to avoid the time costs associated with std::map:
 *     - map<int, list<Vertex>> vertices_for_component_
 *     - map<Vertex,int> component_for_vertex_;
 *     - map<Vertex, vertex_info_t> info_for_vertex;
 *     We're spending 30% time in map::operator[], and probably more in allocation/deletion.
 *     A Vertex is really an integer.
 *
 * S2. dynamic_graph::remove_edge takes 78% of the time.
 *
 * Problems:
 *
 * P1. We don't actually have an implementation of the "fully dynamic graph connectivity"
 *     problem that can quickly delete an edge and determine if this creates a new component!
 *     The HDT approach leads to O(log^2(n)) updates and O(log n/log log) queries.
 *     See: https://people.csail.mit.edu/rahul/papers/dyncon-jea2001.pdf
 *          http://web.stanford.edu/class/archive/cs/cs166/cs166.1166/lectures/17/Small17.pdf
 *
 *     We actually only delete edges, which is the "decremental connectivity" problem, not the "fully dynamic" problem.
 *     Supposedly (1981) https://doi.org/10.1145/322234.322235 implements an O(q + V*E) algorithm.
 *     This approach leads to O(V*E/q) updates and O(1) queries?
 *     See dynamic_graph::remove_edge( ) for an algorithm implementing "process A".
 *     ... except that we are walking the WHOLE computing to remove the smaller piece in lines 415-421.  Not good.
 *     ... also component_for_vertex
 *
 *     Frederickson (1985) https://pdfs.semanticscholar.org/57e3/21514281882ecfa45cc6fc73645e8eecdbbe.pdf
 *     Implements a procedure with O(sqrt(e)) updates and O(1) queries.
 *
 * P2. Problem 1: perf record --call-graph=dwarf is giving an error:
 *     "failed to process sample... failed to process type: 68"
 *     when I run `perf report`.
 *
 * See https://lkml.org/lkml/2016/10/8/189
 *     https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=906728
 */



// bidirectionalS indicates a DIRECTED graph where we can iterate over both in_edges and out_edges
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS> Graph; 
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;

using namespace otc;
namespace fs = boost::filesystem;

using std::vector;
using std::unique_ptr;
using std::set;
using std::list;
using std::tuple;
using std::map;
using std::pair;
using std::string;
using std::optional;
using namespace otc;

// FIXME: I think we want to STOP constructing desids when we read the tree.
//        Maybe we want to record a mapping from each tree node to an integer index...
typedef TreeMappedWithSplits Tree_t;
typedef Tree_t::node_type node_t;

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

struct RSplit {
    vector<int> in;
    vector<int> out;
    vector<int> all;
    RSplit() = default;
    RSplit(const set<int>& i, const set<int>& a) {
        in  = set_to_vector(i);
        all = set_to_vector(a);
        set_difference(begin(all), end(all), begin(in), end(in), std::inserter(out, out.end()));
        assert(in.size() + out.size() == all.size());
    }
};

RSplit split_from_include_exclude(const set<int>& i, const set<int>& e) {
    RSplit s;
    s.in = set_to_vector(i);
    s.out = set_to_vector(e);
    set_union(begin(i),end(i),begin(e),end(e),std::inserter(s.all,s.all.end()));
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
    
    options_description other("Other options");
    other.add_options()
        ("synthesize-taxonomy,T","Make unresolved taxonomy from input tips.")
        ("allow-no-ids,a", "Allow problems w/o OTT ids")
        ("standardize,S", "Write out a standardized subproblem and exit.")
        ;

    options_description visible;
    visible.add(output).add(other).add(otc::standard_options());

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
    write_separated_collection(o, s.in, " ") <<" | ";
    if (s.out.size() < 100) {
        write_separated_collection(o, s.out, " ");
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

bool empty_intersection(const set<int>& xs, const vector<int>& ys)
{
    for (int y: ys){
        if (xs.count(y)) {
            return false;
        }
    }
    return true;
}

struct connected_component_t;

enum class tree_dir {left, right};

tree_dir inverse(tree_dir dir)
{
    if (dir == tree_dir::left)
        return tree_dir::right;
    else
        return tree_dir::left;
}

constexpr unsigned int min_priority = 0;
constexpr unsigned int max_priority = std::numeric_limits<int>::max();

template <typename V>
struct treap_node
{
    static std::mt19937 gen;
    static std::uniform_int_distribution<> priority_dist;

    treap_node<V>* parent = nullptr;
    treap_node<V>* left = nullptr;
    treap_node<V>* right = nullptr;
    const unsigned int priority;
    int n_subtree_nodes = 1;

    const V value;

    treap_node(const V& v)
        :priority(priority_dist(gen)), value(v)
        { }

    treap_node(const V& v, unsigned int p)
        :priority(p), value(v)
        { }
};

template <typename V> std::mt19937 treap_node<V>::gen(0);
template <typename V> std::uniform_int_distribution<> treap_node<V>::priority_dist(min_priority+1,max_priority-1);

// https://en.wikipedia.org/wiki/AA_tree
template <typename V>
struct treap_forest
{
    typedef treap_node<V>* node_t;

    typedef const treap_node<V>* const_node_t;


    bool is_leaf_node(const_node_t node)
    {
        return ((not node->left) and (not node->right));
    }

    bool node_is_valid(const_node_t node)
    {
        if (node->parent and node->parent->priority < node->priority) return false;

        int n = 0;
        if (node->left) n += node->left->n_subtree_nodes;
        if (node->right) n += node->right->n_subtree_nodes;
        if (node->n_subtree_nodes != 1+n) return false;

        return true;
    }

    node_t& get_child(node_t node, tree_dir dir)
    {
        if (dir == tree_dir::left)
            return node->left;
        else
            return node->right;
    }

    void unlink(node_t parent, node_t child, tree_dir dir)
    {
        if (not parent) return;

        assert(parent != child);

        get_child(parent,dir) = nullptr;
        if (child)
        {
            assert(child->parent == parent);
            child->parent = nullptr;
        }
    }

    optional<tree_dir> unlink_from_parent(node_t child)
    {
        auto parent = child->parent;

        if (not parent) return {};

        auto dir = parent_child_dir(parent, child);
        unlink(parent, child, dir);
        return dir;
    }

    // Check that we aren't overwriting any links.
    void link(node_t parent, node_t child, tree_dir dir)
    {
        if (parent)
        {
            assert(parent != child);
            // there is not already a child here
            assert(not get_child(parent,dir));
            get_child(parent,dir) = child;
            if (child)
                child->parent = parent;
        }
        else if (child)
        {
            assert(not child->parent);
        }
    }

    void swap_subtrees(node_t sub1, node_t sub2)
    {
        assert(sub2);
        assert(not sub2->parent);

        auto parent1 = sub1->parent;
        auto dir1 = unlink_from_parent(sub1);

        auto parent2 = sub2->parent;
        auto dir2 = unlink_from_parent(sub2);

        if (dir1)
            link(parent1, sub2, *dir1);

        if (dir2)
            link(parent2, sub2, *dir2);
    }

    void update_subtree_nodes_one(node_t node)
    {
        auto& n = node->n_subtree_nodes;
        n = 1;
        if (node->left) n += node->left->n_subtree_nodes;
        if (node->right) n += node->right->n_subtree_nodes;
    }

    void update_subtree_nodes(node_t node)
    {
        while (node)
        {
            update_subtree_nodes_one(node);
            node = node->parent;
        }
    }

    tree_dir parent_child_dir(node_t parent, node_t child)
    {
        assert(parent->left == child or parent->right == child);
        if (parent->left == child)
            return tree_dir::left;
        else if (parent->right == child)
            return tree_dir::right;
        std::abort();
    }

    /*
    //   grandparent           grandparent
    //       /                      /
    //    parent                 child
    //    //  \          ->      /   \\
    //    X   child           parent  Y
    //        /   \\          //  \
    //       A     Y          X    A
    */

    // Moves the child up to the parent's position
    void rotate(node_t parent, node_t child)
    {
        assert(child);
        assert(child->parent);
        assert(child->parent == parent);

        // grandparent could be NULL.
        auto grandparent = parent->parent;
        tree_dir parent_dir;
        if (grandparent)
            parent_dir = parent_child_dir(grandparent, parent);
        auto child_dir = parent_child_dir(parent, child);
        auto A_dir = inverse(child_dir);
        node_t A = get_child(child, A_dir);

#ifndef NDEBUG
        auto x = parent->n_subtree_nodes;
#endif

        unlink(grandparent, parent, parent_dir);
        unlink(parent, child, child_dir);
        unlink(child, A, A_dir);

        link(parent, A, child_dir);
        link(child, parent, A_dir);
        link(grandparent, child, parent_dir);

        update_subtree_nodes_one(parent);
        update_subtree_nodes_one(child);

        assert(x == child->n_subtree_nodes);
    }

    node_t last_in_subtree(node_t v1)
    {
        assert(v1);
        while(v1->right)
            v1 = v1->right;
        return v1;
    }

    node_t prev_from_subtree(node_t node)
    {
        assert(node);
        if (node->left)
            return last_in_subtree(node->left);
        else
            return nullptr;
    }

    node_t first_in_subtree(node_t v1)
    {
        assert(v1);
        while(v1->left)
            v1 = v1->left;
        return v1;
    }

    node_t next_from_subtree(node_t node)
    {
        assert(node);
        if (node->right)
            return first_in_subtree(node->right);
        else
            return nullptr;
    }

    // Find the first node in a tour
    node_t first_in_tour(node_t v1)
    {
        return first_in_subtree(root(v1));
    }

    // Find the last node in a tour
    node_t last_in_tour(node_t v1)
    {
        return last_in_subtree(root(v1));
    }

    node_t prev(node_t node)
    {
        if (auto p = prev_from_subtree(node))
            return p;

        while (auto parent = node->parent)
        {
            if (parent->right == node)
                return parent;

            node = parent;
        }

        return nullptr;
    }

    node_t next(node_t node)
    {
        if (auto n = next_from_subtree(node))
            return n;

        while (auto parent = node->parent)
        {
            if (parent->left == node)
                return parent;

            node = parent;
        }

        return nullptr;
    }

    node_t next(node_t node, tree_dir dir)
    {
        if (dir == tree_dir::right)
            return next(node);
        else
            return prev(node);
    }

    // Find the root node in a tour
    node_t root(node_t v1)
    {
        assert(v1);
        while (v1->parent)
            v1 = v1->parent;
        return v1;
    }

    int size_of_component(node_t node)
    {
        if (node)
            return root(node)->n_subtree_nodes;
        else
            return 0;
    }

    bool same_component(node_t node1, node_t node2)
    {
        assert(node1);
        assert(node2);
        return root(node1) == root(node2);
    }

    std::uint64_t directions_to_root(node_t node)
    {
        std::uint64_t path = 0;
        std::uint64_t bit = 1UL<<63;
        while(auto parent = node->parent)
        {
            if (parent->right == node)
                path |= bit;
            path >>= 1;
            node = parent;
        }
        return path;
    }

    bool is_smaller(node_t u, node_t v)
    {
        assert(root(u) == root(v));

        return directions_to_root(u) < directions_to_root(v);
    }

private:
    void place_node(node_t pos, tree_dir dir, node_t node)
    {
        assert(pos);
        assert(node);
        // 2. Insert the node at the correct position
        if (auto pos2 = get_child(pos,dir))
        {
            if (dir == tree_dir::left)
                pos2 = last_in_subtree(pos2);
            else
                pos2 = first_in_subtree(pos2);
            link(pos2, node, inverse(dir));
        }
        else
            link(pos, node, dir);
    }

public:
    void insert(node_t pos, tree_dir dir, node_t node)
    {
        assert(node);

        // 1. If the treap is empty, just return the node.
        if (not pos) return;

        // 2. Add the node
        place_node(pos, dir, node);

        // 3. Update the number of subtree nodes in all ancestors
        update_subtree_nodes(node->parent);

        // 4. Rebalance the treap
        while(node->parent and node->parent->priority < node->priority)
            rotate(node->parent, node);
    }

    node_t insert(node_t pos, tree_dir dir, const V& value)
    {
        auto node = new treap_node<V>(value);
        insert(pos, dir, node);
        return node;
    }

    node_t isolate(node_t node)
    {
        assert(node);

        // 1. Rotate node down until it has < 2 children.
        while (node->left and node->right)
        {
            auto child = node->left;
            if (node->right->priority > node->left->priority)
                child = node->right;
            rotate(node,child);
        }

        // 2. Unlink the node if it is not a leaf.
        auto parent = node->parent;
        node_t other = parent;
        if (auto left = node->left)
        {
            if (not parent) other = left;
            unlink(node, left, tree_dir::left);
            swap_subtrees(node, left);
        }
        else if (auto right = node->right)
        {
            if (not parent) other = right;
            unlink(node, right, tree_dir::right);
            swap_subtrees(node, right);
        }
        else
            unlink_from_parent(node);

        // 3. Update subtree nodes count
        update_subtree_nodes(parent);
        return other;
    }

    node_t remove(node_t node)
    {
        auto other = isolate(node);
        delete node;
        return other;
    }

    pair<node_t,node_t> split(node_t node, tree_dir dir)
    {
        // 1. Make a dummy tree node
        treap_node<V> _dummy({}, max_priority);
        node_t dummy = &_dummy;

        insert(node, dir, dummy);
        auto left = dummy->left;
        auto right = dummy->right;

        unlink(dummy, left , tree_dir::left );
        unlink(dummy, right, tree_dir::right);

        return {left, right};
    }

    node_t join(node_t left, node_t right)
    {
        // 1. If one treap is empty, just return the other.
        if (not right)
            return left;
        if (not left)
            return right;

        // 2. Find the roots of the two trees
        left = root(left);

        right = root(right);

        assert(left != right);

        // 3. Make a dummy tree node
        treap_node<V> _dummy({}, min_priority);
        node_t dummy = &_dummy;

        // 4. Connect dummy to left and right subtrees
        link(dummy, left, tree_dir::left);
        link(dummy, right, tree_dir::right);
        update_subtree_nodes_one(dummy);

        // 5. Remove the dummy from the treap
        isolate(dummy);

        return left;
    }

    node_t make_first(node_t v1)
    {
        if (not v1) return nullptr;

        auto r = root(v1);
        auto w1 = first_in_tour(r);

        // If v1 is already first, then quit here
        if (v1 == w1) return r;

        auto [prefix,postfix] = split(v1, tree_dir::left);

        return join(postfix,prefix);
    }

    string show_treap(node_t v);
};

template <typename V>
std::ostream& operator<<(std::ostream& o, treap_node<V>* v)
{
    treap_forest<V> F;
    o<<" ";
    v = F.first_in_tour(v);
    while(v)
    {
        assert(F.node_is_valid(v));
        o<<v->value<<" ";
        v = F.next(v);
    }
    return o;
}

template <typename V>
string treap_forest<V>::show_treap(treap_forest<V>::node_t node)
{
    std::ostringstream o;
    o<<node;
    return o.str();
}

class vertex_info_t
{
    // We store a pointer to the component for every marked node.
    //   To test if a node is marked, we test if the pointer is non-null.
    // Nodes are marked if the have no parent.  Such nodes are precisely the
    //   nodes in the in a connected component Y (of the display graph)
    //   that make up the corresponding "position" U (of the cluster graph).
    connected_component_t* component;

public:
    // Which tree (0..k-1) is this vertex in?
    // If the vertex corresponds to a label, then it is not in a tree, and the index is unset.
    // (This differs from the paper, which uses 1..k for trees, and k+1 for no tree.)
    optional<int> tree_index;

    optional<OttId> label;

    mutable char flags;

    int component_index;

    std::list<Vertex>::iterator list_entry;

    bool is_marked() const {return component;}
    void unmark_node() {component = nullptr;}
    void mark_node(connected_component_t* c) {component = c;}

    connected_component_t* get_component() const {assert(component); return component;}
};

/*
 * Q: Do we need to remove the map from each vertex to its component index?
 * A: ??
 *    Probably, since then splitting a component would require walking every vertex in the smaller piece?
 *
 * Q: How do we know if each edge is a tree edge or a non-tree edge?
 * A: ??
 *    We could put all tree edges in a tree_edge hash?
 *
 * Q: Can we store component info on the root of each spanning tree?
 *
 * 1. OK, so every vertex is part of a spanning tree.
 * 2. We need to know if each edge is a tree edge or not
 *
 */
template<typename T> void hash_combine(size_t & seed, T const& v)
{
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct edge
{
    long unsigned int v1;
    long unsigned int v2;
public:
    long unsigned int source() const {return v1;}
    long unsigned int target() const {return v2;}

    bool operator==(const edge& e2) const { return v1 == e2.v1 and v2 == e2.v2;}

    [[nodiscard]] edge reverse() const;

    edge(long unsigned int x1, long unsigned int x2)
        :v1(x1),v2(x2)
        { }

    edge():edge(-1,-1) {}
};

std::ostream& operator<<(std::ostream& o, const edge& e)
{
    o<<"{"<<e.source()<<","<<e.target()<<"}";
    return o;
}

edge edge::reverse() const
{
    return edge(target(), source());
}

namespace std
{
    template<>
    class hash<edge> {
    public:
        size_t operator()(const edge &e) const
            {
                size_t seed = std::hash<int>()(e.source());
                hash_combine(seed, std::hash<int>()(e.target()));
                return seed;
            }
    };
}


enum edge_type_t {tree_edge = 0, non_tree_edge = 1};

struct edge_info_t
{
    edge_type_t type;
    treap_node<edge>* euler_tour_node;
};


class dynamic_graph
{
    Graph G;

    map<int,list<Vertex>> vertices_for_component_;

    int next_component = 0;

    vector<vertex_info_t> info_for_vertex;

    // Probably this should be an edge property.
    robin_hood::unordered_map<edge, edge_info_t> edge_info;

    static treap_forest<edge> F;

    typedef treap_node<edge>* euler_tour_node_t;

    typedef const treap_node<edge>* const_euler_tour_node_t;

public:
    const vertex_info_t& vertex_info(Vertex v) const {return info_for_vertex[v];}

          vertex_info_t& vertex_info(Vertex v)       {return info_for_vertex[v];}

    int new_component() {return next_component++;}

    char& flags(Vertex v) const {return vertex_info(v).flags;}

    int n_components() const {return vertices_for_component_.size();}

    int& component_for_vertex(Vertex v)       {return vertex_info(v).component_index;}

    int  component_for_vertex(Vertex v) const {return vertex_info(v).component_index;}

    const map<int,list<Vertex>>& vertices_for_components() const {return vertices_for_component_;}

          map<int,list<Vertex>>& vertices_for_components()       {return vertices_for_component_;}

    const list<Vertex>& vertices_for_component(int c) const {return vertices_for_component_.at(c);}

          list<Vertex>& vertices_for_component(int c)       {return vertices_for_component_.at(c);}

    int size_of_component(int c) const {return vertices_for_component(c).size();}

    auto in_edges(Vertex v) const {return boost::in_edges(v,G);}

    auto out_edges(Vertex v) const {return boost::out_edges(v,G);}

    int num_vertices() const {return boost::num_vertices(G);}

    int num_edges() const {return boost::num_edges(G);}

    Vertex source(Edge e) const {return boost::source(e,G);}

    Vertex target(Edge e) const {return boost::target(e,G);}

    int in_degree(Vertex v) const
    {
        int count = 0;
        for(auto [e, e_end] = in_edges(v); e != e_end; e++)
            count++;
        return count;
    }

    int out_degree(Vertex v) const
    {
        int count = 0;
        for(auto [e, e_end] = out_edges(v); e != e_end; e++)
            count++;
        return count;
    }

    bool is_tip_node(Vertex v) const;

    bool is_internal_node(Vertex v) const
    {
        return not is_tip_node(v);
    }

    Vertex add_vertex();

    bool same_component(Vertex u, Vertex v) const
    {
        int cu = component_for_vertex(u);
        int cv = component_for_vertex(v);

        return (cu == cv);
    }

    euler_tour_node_t to_euler_tour_node(const optional<Edge> E) const
    {
        if (E)
        {
            edge e{source(*E),target(*E)};
            return edge_info.at(e).euler_tour_node;
        }
        else
            return nullptr;
    }

    euler_tour_node_t to_reverse_euler_tour_node(const optional<Edge> E) const
    {
        if (E)
        {
            edge e{target(*E),source(*E)};
            return edge_info.at(e).euler_tour_node;
        }
        else
            return nullptr;
    }

    optional<Edge> some_tree_edge_from(Vertex u) const
    {
        for(auto [e, e_end] = out_edges(u); e != e_end; e++)
            if (is_tree_edge(*e))
                return *e;
        return {};
    }

    optional<Edge> some_tree_edge_to(Vertex u) const
    {
        for(auto [e, e_end] = in_edges(u); e != e_end; e++)
            if (is_tree_edge(*e))
                return *e;
        return {};
    }

    euler_tour_node_t some_node_to(Vertex v) const
    {
        if (auto E = some_tree_edge_to(v))
            return to_euler_tour_node(*E);
        else
            return to_reverse_euler_tour_node(some_tree_edge_from(v));
    }

    euler_tour_node_t some_node_from(Vertex v) const
    {
        if (auto E = some_tree_edge_from(v))
            return to_euler_tour_node(*E);
        else
            return to_reverse_euler_tour_node(some_tree_edge_to(v));
    }
    
    bool same_spanning_tree(Vertex u, Vertex v) const
    {
        // 1. If u and v are the same vertex, then return true.
        if (u == v) return true;

        // 2. Get pointers to Euler tour entries for each component
        auto node_u = some_node_to(u);

        auto node_v = some_node_to(v);

        // 3. If vertex has no Euler tour entries, then that vertex has no edges, and is in its own unique component.
        if (not node_u or not node_v)
            return false;

        // 4. If both vertices have adjacent edges, then check if they are in the same Euler tour tree.
        return F.same_component(node_u, node_v);
    }

    string show_tour(euler_tour_node_t node) const;

    bool component_smaller(Vertex u, Vertex v) const
    {
        int cu = component_for_vertex(u);
        int cv = component_for_vertex(v);

        return size_of_component(cu) < size_of_component(cv);
    }

    int size_of_component2(Vertex u) const
    {
        auto e = some_node_from(u);
        auto n_edges = F.size_of_component(e);
        assert(n_edges%2 == 0);
        return 1 + (n_edges/2);
    }

    bool component_smaller2(Vertex u, Vertex v) const
    {
        return size_of_component2(u) < size_of_component2(v);
    }

    void add_tree_edge(Vertex u, Vertex v)
    {
        edge E1(u,v);
        edge E2(v,u);
        assert(not edge_info.count(E1) or edge_info.at(E1).type == edge_type_t::non_tree_edge);
        assert(not edge_info.count(E2) or edge_info.at(E2).type == edge_type_t::non_tree_edge);

        // Finding Eu and Ev has to happen before we add the (u,v) edge to the graph.
        auto Eu = some_node_from(u);
        F.make_first(Eu);

        auto Ev = some_node_from(v);
        F.make_first(Ev);

        auto node_uv = F.insert(nullptr, tree_dir::right, E1);
        auto node_vu = F.insert(nullptr, tree_dir::right, E2);
        // This overwrites any previous status.  Previously it might have been a non-tree edge.
        edge_info[E1] = edge_info_t{tree_edge, node_uv};
        edge_info[E2] = edge_info_t{tree_edge, node_vu};

        auto Euv = F.join(F.join(F.join(Eu,node_uv),Ev), node_vu);
    }

    auto add_edge(Vertex u, Vertex v)
    {
        // 1. Check that we don't already have an edge (u,v) or (v,u)
        {
            auto [_,found] = boost::edge(u,v,G);
            assert(not found);
            auto [_2,found2] = boost::edge(v,u,G);
            assert(not found2);
        }

        // 2. If the vertices are in different components, it is a tree edge.
        //    Otherwise, it is a non-tree edge.
        bool same_comp = same_component(u,v);
        bool same_comp2 = same_spanning_tree(u,v);
        assert(same_comp == same_comp2);

        assert(size_of_component(component_for_vertex(u)) == size_of_component2(u));
        assert(size_of_component(component_for_vertex(v)) == size_of_component2(v));
        assert(component_smaller(v,u) == component_smaller2(v,u));

        // 3a. Add a non-tree edge if both vertices in same component
        if (same_comp)
        {
            edge E1(u,v);
            edge E2(v,u);
            assert(not edge_info.count(E1));
            assert(not edge_info.count(E2));

            edge_info.insert({E1, edge_info_t{non_tree_edge,nullptr}});
            edge_info.insert({E2, edge_info_t{non_tree_edge,nullptr}});
            auto e = boost::add_edge(u,v,G);
            return e;
        }
        // 3b. Merge the components if they are different
        else
        {
            // Join the two euler tour paths.
            add_tree_edge(u,v);

            auto e = boost::add_edge(u,v,G);

            // 2a. Ensure that cu is smaller, or equal.
            if (component_smaller(v,u)) std::swap(u,v);

            int cu = component_for_vertex(u);
            int cv = component_for_vertex(v);

            int newsize = size_of_component(cu)+size_of_component(cv);

            auto& lu = vertices_for_component(cu);
            auto& lv = vertices_for_component(cv);

            // 2b. Mark all the vertices in cu as no being in cv
            for(auto uu: lu)
                component_for_vertex(uu) = cv;

            // 2c. Move the vertex list from cu to the cv component.
            lv.splice(lv.end(), lu);

            // 2d. There is no longer a cu component.
            vertices_for_component_.erase(cu);

            assert(same_component(u,v));
            assert(same_spanning_tree(u,v));

            assert(size_of_component2(u) == newsize);
            assert(size_of_component2(v) == newsize);
            assert(size_of_component(cv) == newsize);
            return e;
        }
    }

    bool is_tree_edge(const edge& e) const
    {
        edge E(e.source(), e.target());
        assert(edge_info.count(E));
        return edge_info.at(E).type == tree_edge;
    }

    bool is_tree_edge(Vertex u, Vertex v) const
    {
        return is_tree_edge(edge(u,v));
    }

    bool is_tree_edge(const Edge& e) const
    {
        return is_tree_edge(source(e), target(e));
    }

    vector<Vertex> find_spanning_tree_for_vertex2(Vertex u) const
    {
        auto node = some_node_from(u);
        if (not node)
            return {u};

        vector<Vertex> vertices;

        // For each edge in the tour, consider its target node.
        for(node = F.first_in_tour(node); node; node = F.next(node))
        {
            auto v = node->value.target();
            if (not flags(v))
            {
                flags(v) = 1;
                vertices.push_back(v);
            }
        }

        for(auto& v: vertices)
            flags(v) = 0;

        return vertices;
    }

    vector<Vertex> find_spanning_tree_for_vertex(Vertex u) const
    {
        vector<Vertex> nodes;

        // Return true if vertices are in the same component
        auto try_add = [&](Vertex uu)
                           {
                               if (not flags(uu))
                               {
                                   flags(uu) = 1;
                                   nodes.push_back(uu);
                               }
                           };

        try_add(u);

        int i=0;
        for(;i<nodes.size();i++)
        {
            auto uu = nodes[i];

            for(auto [e, e_end] = in_edges(uu); e != e_end; e++)
                if (is_tree_edge(*e))
                    try_add( source(*e) );

            for(auto [e, e_end] = out_edges(uu); e != e_end; e++)
                if (is_tree_edge(*e))
                    try_add( target(*e) );
        }

        for(auto& node: nodes)
            flags(node) = 0;

        return nodes;
    }

    pair<euler_tour_node_t, euler_tour_node_t>
    remove_tree_edge(Vertex u, Vertex v)
    {
        // 1. Find the nodes for (u,v) and (v,u) in the tour.
        edge edge_uv(u,v);
        edge edge_vu(v,u);
        auto node_uv = edge_info.at(edge_uv).euler_tour_node;
        auto node_vu = edge_info.at(edge_vu).euler_tour_node;
        assert(node_uv);
        assert(node_vu);

        // 2. Make the tour start at (u,v)
        F.make_first(node_uv);

        // 3. Split into a tour beginning with (u,v) and a tour beginning with (v,u)
        auto [Ev,Eu] =  F.split(node_vu, tree_dir::left);

        // 4a. Remove (u,v) to get the tour beginning and ending with v
        Ev = F.remove(node_uv);
        // 4b. Remove (v,u) to get the tour beginning and ending with u
        Eu = F.remove(node_vu);

        // 5. Remove edge annotations.
        edge_info.erase(edge_uv);
        edge_info.erase(edge_vu);

        return {Eu,Ev};
    }

    // if the edge is a non-tree edge:
    //    we don't have to do anything
    // else
    //    cut the tree edge
    //    find a replacement non-tree edge if one exists
    //    mark it as a tree edge
    bool remove_edge(Vertex u, Vertex v)
    {
        edge E(u,v);
        bool was_tree_edge = is_tree_edge(E);

        assert(same_spanning_tree(u,v));

        if (was_tree_edge)
            remove_tree_edge(u,v);

        for(int v = 0; v < num_vertices(); v++)
            assert(flags(v) == 0);

        int c1 = component_for_vertex(u);
        boost::remove_edge(u,v,G);
        vector<Vertex> from_u;
        vector<Vertex> from_v;

        // Return true if vertices are in the same component
        auto try_add_u = [&](Vertex uu)
                             {
                                 if (not flags(uu))
                                 {
                                     flags(uu) = 1;
                                     from_u.push_back(uu);
                                     return false;
                                 }
                                 else if (flags(uu) != 1)
                                     return true;
                                 else
                                     return false;
                             };

        auto try_add_v = [&](Vertex vv)
                             {
                                 if (not flags(vv))
                                 {
                                     flags(vv) = 2;
                                     from_v.push_back(vv);
                                     return false;
                                 }
                                 else if (flags(vv) != 2)
                                     return true;
                                 else
                                     return false;
                             };

        try_add_u(u);
        try_add_v(v);

        int i=0,j=0;

        bool same_component = false;

        while(i < from_u.size() and j < from_v.size())
        {
            // Check 1 entry from u
            auto uu = from_u[i++];
            for(auto [e, e_end] = in_edges(uu); e != e_end; e++)
                if (try_add_u( source( *e ) ))
                {
                    same_component = true;
                    break;
                }
            for(auto [e, e_end] = out_edges(uu); e != e_end; e++)
                if (try_add_u( target( *e ) ))
                {
                    same_component = true;
                    break;
                }

            // Check 1 entry from v
            auto vv = from_v[j++];
            for(auto [e, e_end] = in_edges(vv); e != e_end; e++)
                if (try_add_v( source( *e ) ))
                {
                    same_component = true;
                    break;
                }
            for(auto [e, e_end] = out_edges(vv); e != e_end; e++)
                if (try_add_v( target( *e ) ))
                {
                    same_component = true;
                    break;
                }
        }

        // Clear the seen flags
        for(auto u: from_u)
            flags(u) = 0;

        for(auto v: from_v)
            flags(v) = 0;

        for(int v = 0; v < num_vertices(); v++)
            assert(flags(v) == 0);

        optional<edge> connecting_tree_edge;
        if (was_tree_edge)
        {
            auto spanning_tree_for_u = find_spanning_tree_for_vertex(u);
            auto spanning_tree_for_u2 = find_spanning_tree_for_vertex2(u);
            assert(spanning_tree_for_u.size() == spanning_tree_for_u2.size());
            auto spanning_tree_for_v = find_spanning_tree_for_vertex(v);
            assert(spanning_tree_for_u.size() + spanning_tree_for_v.size() == vertices_for_component(c1).size());

            if (spanning_tree_for_u.size() > spanning_tree_for_v.size())
            {
                std::swap(u,v);
                std::swap(spanning_tree_for_u, spanning_tree_for_v);
            }

            for(auto& uu: spanning_tree_for_u)
                assert(flags(uu) == 0);
            for(auto& uu: spanning_tree_for_u)
                flags(uu) = 1;
            for(auto& vv: spanning_tree_for_v)
                assert(flags(vv) == 0);

            for(auto& uu: spanning_tree_for_u)
            {
                for(auto [e, e_end] = in_edges(uu); e != e_end; e++)
                {
                    if (not is_tree_edge(*e) and flags(source(*e)) == 0)
                    {
                        connecting_tree_edge = edge(source(*e), target(*e));
                        break;
                    }
                }
                for(auto [e, e_end] = out_edges(uu); e != e_end; e++)
                {
                    if (not is_tree_edge(*e) and flags(target(*e)) == 0)
                    {
                        connecting_tree_edge = edge(source(*e), target(*e));
                        break;
                    }
                }
            }

            for(auto& uu: spanning_tree_for_u)
                flags(uu) = 0;
        }

        // Quit here if we didn't split a component
        if (same_component)
        {
            assert(not was_tree_edge or connecting_tree_edge);
            if (was_tree_edge)
            {
                assert(connecting_tree_edge);
                auto edge_wx = *connecting_tree_edge;

                assert(edge_info.count(edge_wx));
                assert(edge_info.count(edge_wx.reverse()));

                auto w = edge_wx.source();
                auto x = edge_wx.target();

                add_tree_edge(w,x);
            }
            return false;
        }
        else
        {
            assert(was_tree_edge and not connecting_tree_edge);
        }

        // Move vertices from the smaller group to a new component
        if (i < from_u.size())
            std::swap(from_u, from_v);

        int c2 = new_component();
        auto& nodes1 = vertices_for_component(c1);
        list<Vertex> nodes2;
        for(auto uu: from_u)
        {
            component_for_vertex(uu) = c2;
            auto it = vertex_info(uu).list_entry;
            nodes2.splice(nodes2.end(), nodes1, it);
        }
        vertices_for_component_[c2] = std::move(nodes2);

        return true;
    }
};

treap_forest<edge> dynamic_graph::F;

string dynamic_graph::show_tour(euler_tour_node_t node) const
{
    return F.show_treap(node);
}

bool semi_universal_position(const dynamic_graph& G, const set<Vertex>& vs)
{
    if (vs.size() == 1)
    {
        auto v = *vs.begin();
        if (G.is_internal_node(v)) return true;
    }
    return false;
}

bool dynamic_graph::is_tip_node(Vertex v) const
{
    assert(vertex_info(v).tree_index);
    auto [e, e_end] = out_edges(v);
    
    // We shouldn't have any out-degree-0 nodes.
    assert(e != e_end);
    auto f = e;
    f++;
    if (f == e_end)
    {
        // Degree-1 node
        assert(vertex_info(boost::target(*e,G)).label);
        return true;
    }
    else
        return false;
}

Vertex dynamic_graph::add_vertex()
{
    auto v = boost::add_vertex(G);
    info_for_vertex.push_back({});
    int c = new_component();
    component_for_vertex(v) = c;
    vertices_for_component_[c] = {v};
    assert(component_for_vertex(v) == c);

    vertex_info(v).list_entry = vertices_for_component_[c].begin();

    assert(*vertex_info(v).list_entry == v);
    assert(info_for_vertex.size() == num_vertices());

    return v;
}

struct connected_component_t
{
    dynamic_graph* G = nullptr;

    dynamic_graph* get_graph() {return G;}

    // Y.count: the cardinality of Y \cap X_P
    // That is, the number of tip labels in the component.
    int count = 0;

    // Y.map: {(i,L[i])| i is a tree index, and L[i] is not empty}
    //        where L[i] is the set of marked vertices in tree T[i].
    map<int,set<Vertex>> marked_vertices_for_tree;

    void erase_marked_vertex(Vertex v)
    {
        assert(G->vertex_info(v).tree_index);
        int i = *G->vertex_info(v).tree_index;
        assert(G->vertex_info(v).is_marked());
        auto& nodes = marked_vertices_for_tree.at(i);
        if (semi_universal_position(*G, nodes))
            semiU.erase(i);
        nodes.erase(v);
        if (semi_universal_position(*G, nodes))
            semiU.insert(i);
        G->vertex_info(v).unmark_node();

        // If there are no marked vertices left, then erase the record.
        if (nodes.size() == 0)
            marked_vertices_for_tree.erase(i);
    }

    void insert_marked_vertex(Vertex v)
    {
        assert(G->vertex_info(v).tree_index);
        int i = *G->vertex_info(v).tree_index;
        assert(not G->vertex_info(v).is_marked());
        auto& nodes = marked_vertices_for_tree[i];
        if (semi_universal_position(*G, nodes))
            semiU.erase(i);
        nodes.insert(v);
        if (semi_universal_position(*G, nodes))
            semiU.insert(i);
        G->vertex_info(v).mark_node(this);
    }

    // Y.semiU: the indices i for trees where Y.map[i] is defined
    //          (i.e. L[i] is not empty) AND |Y.map[i]| == 1, AND
    //          the single element of Y.map[i] is an internal node in T[i].
    set<int> semiU;

    vector<Vertex> semi_universal_nodes_for_position() const
    {
        vector<Vertex> SU;
        for(int index: semiU)
            for(Vertex v: marked_vertices_for_tree.at(index))
                SU.push_back(v);
        return SU;
    }

    vector<Vertex> vertices_in_component() const
    {
        vector<Vertex> vertices;
        set<Vertex> visited;

        // Add initial vertices
        for(auto& [index,vs]: marked_vertices_for_tree)
            for(auto& v: vs)
            {
                vertices.push_back(v);
                visited.insert(v);
            }

        for(int i=0;i<vertices.size();i++)
        {
            for(auto [e, e_end] = G->in_edges(vertices[i]); e != e_end; e++)
            {
                auto v = G->source( *e );
                if (not visited.count(v))
                {
                    visited.insert(v);
                    vertices.push_back(v);
                }
            }
            for(auto [e, e_end] = G->out_edges(vertices[i]); e != e_end; e++)
            {
                auto v = G->target( *e );
                if (not visited.count(v))
                {
                    visited.insert(v);
                    vertices.push_back(v);
                }
            }
        }

        return vertices;
    }

    vector<OttId> labels_for_component() const
    {
        vector<OttId> labels;
        auto vs = vertices_in_component();
        for(auto v: vs)
            if (auto l = G->vertex_info(v).label)
                labels.push_back(*l);
        return labels;
    }

    connected_component_t(dynamic_graph* g):G(g) {}
};

typedef vector<const node_t*> position_t;

void append_leaves_not_under(const node_t* start, const node_t* avoid, vector<const node_t*>& leaves)
{
    // 1. Don't consider the node we are avoiding
    if (start == avoid) return;

    // 2. If this is a tip, then count it.
    if (start->is_tip())
    {
        leaves.push_back(start);
        return;
    }

    // 3. If this is NOT a tip, consider all of its children
    for(auto node = start->get_first_child(); node ; node = node->get_next_sib())
        append_leaves_not_under(node, avoid, leaves);
}

vector<const node_t*> leaves_not_under(const node_t* avoid)
{
    // 1. Find the true root.
    auto root = avoid;
    while(root->get_parent())
        root = root->get_parent();

    // 2. Find leaves under the root that are not under avoid.
    vector<const node_t*> leaves;
    append_leaves_not_under(root, avoid, leaves);

    return leaves;
}


pair<unique_ptr<dynamic_graph>,vector<unique_ptr<connected_component_t>>>
display_graph_from_profile(const vector<const node_t*>& profile)
{
    unique_ptr<dynamic_graph> H(new dynamic_graph);
    map<OttId,Vertex> Labels;
    map<const node_t*,Vertex> node_to_vertex;

    auto vertex_for_label = [&](const node_t* nd)
                                {
                                    assert(nd->has_ott_id());
                                    auto id = nd->get_ott_id();

                                    if (not Labels.count(id))
                                    {
                                        auto label = H->add_vertex();
                                        Labels.insert({id,label});
                                        H->vertex_info(label).label = id;
                                    }

                                    return Labels.at(id);
                                };

    // 1. Construct the display graph
    for(int i=0; i< profile.size(); i++)
    {
        // The root should be an actual root now (no parent).
        auto root = profile[i];
        assert(not root->get_parent());

        // FIXME!  How should we handle incertae sedis taxa?

        // Walk nodes in the subtree below profile[i] in pre-order (parent before child) so that we can connect children to parents.
        for(auto nd: iter_pre_n_const(root))
        {
            // Add vertex for child node.
            auto v = H->add_vertex();
            node_to_vertex.insert({nd,v});
            H->vertex_info(v).tree_index = i;

            // Add edge from parent node to child node
            if (nd != root)
            {
                auto u = node_to_vertex.at(nd->get_parent());
                H->add_edge(u,v);
            }

            // Add an edge from the node for the label to the tip
            if (nd->is_tip())
                H->add_edge(v, vertex_for_label(nd));
        }
    }

    // 2. Construct the Y_inits

    // 2a. Function to find which component should contain a particular tree root.
    vector<unique_ptr<connected_component_t>> Y_inits;
    auto find_component = [&](Vertex u) -> optional<int>
    {
        for(int i=0;i<Y_inits.size();i++)
        {
            auto& Y_init = Y_inits[i];

            // 1. Find a vertex in the component
            optional<Vertex> v;
            for(auto& [tree_index,vertices]: Y_init->marked_vertices_for_tree)
            {
                assert(not vertices.empty());
                v = *vertices.begin();
                break;
            }
            assert(v);

            // 2. If u is connected to v, then u is in component i.
            if (H->same_spanning_tree(u, *v))
                return i;
        }
        return {};
    };

    // 2b. Find the component for each tree
    vector<int> component_for_tree;
    for(int i=0; i< profile.size(); i++)
    {
        auto root_vertex = node_to_vertex.at(profile[i]);
        auto component_index = find_component(root_vertex);
        if (not component_index)
        {
            component_index = Y_inits.size();
            Y_inits.push_back(unique_ptr<connected_component_t>(new connected_component_t(&*H)));
        }
        Y_inits[*component_index]->insert_marked_vertex(root_vertex);
        component_for_tree.push_back( *component_index );
    }

    // 2c. Assign each label to a component
    for(auto& [ottid,vertex]: Labels)
    {
        auto c = find_component(vertex);
        assert(c);
        Y_inits[*c]->count++;
    }

    
    LOG(DEBUG)<<"Labels";
    for(auto [id,_]: Labels)
    {
        LOG(DEBUG)<<"   "<<id;
    }
    LOG(DEBUG)<<"";

    return pair{std::move(H),std::move(Y_inits)};
}

// Walk the component Y2 containing node. Move marked nodes from Y1 to Y2.
void split_component(connected_component_t* Y1, connected_component_t* Y2, Vertex node1)
{
    auto G = Y1->get_graph();
    int cnode1 = G->component_for_vertex(node1);

    assert(Y1->count >  0);
    assert(Y2->count == 0);
    
    for(auto node2: G->vertices_for_component(cnode1))
    {
        auto& info = G->vertex_info(node2);
        if (not info.tree_index)
        {
            Y1->count--;
            Y2->count++;
            assert(Y1->count > 0);
            continue;
        }
        int i = *info.tree_index;
        if (info.is_marked())
        {
            Y1->erase_marked_vertex(node2);
            Y2->insert_marked_vertex(node2);
        }
    }
}

// We need to return a new component corresponding to the second vertex u
unique_ptr<connected_component_t> split_component(connected_component_t* Y1, Vertex parent, Vertex child)
{
    auto G = Y1->get_graph();
    assert(not G->vertex_info(parent).is_marked());
    assert(G->vertex_info(child).get_component() == Y1);

    unique_ptr<connected_component_t> W(new connected_component_t(Y1->get_graph()));
    auto Y2 = W.get();

//    assert(Y1->count == Y1->labels_for_component().size());
    int cparent = G->component_for_vertex(parent);
    int cchild = G->component_for_vertex(child);

    if (G->size_of_component(cchild) > G->size_of_component(cparent))
        // Move marked nodes connected to parent (parent) to Y2
        split_component(Y1, Y2, parent);
    else
        // Move marked nodes connected to child (child) to Y2
        split_component(Y1, Y2, child);

//    assert(Y1->count == Y1->labels_for_component().size());
//    assert(Y2->count == Y2->labels_for_component().size());

    // If we are disconnecting a parent from its last child, then Y2
    // could be a component containing only the parent node and no labels.
    //
    // However, in that case, we should be able to keep Y1 as representing
    // the child.
    assert(Y1->count != 0);

    if (W->count > 0)
        return W;
    else
        return {};
}

// There will always initially be one component containing all vertices...?
// The paper says this.
// But what if there are only two trees, and they have non-overlapping leaf sets?

// Question: How do we get from the positions U to the connected components Y?
// Answer:   According to Lemma 4, the connected components of a G_P(U) *are* positions AND valid positions.
//           U is just a set of nodes on the profile, and G_P(u) connects vertices u and v
//           if cluster(u) intersects cluster(v).
//

// Question: How about for H_P(U)?
// Answer:   H_P(U) is the subgraph induced by keeping only descendants of nodes v \in U.
//           Roots of H_P(U) correspond to nodes in U.
//
//           From section 4.4, the connected components W[i] of G_P(U) can be put into a 1-1
//           correspondence with the connected components Y[i] of H_P(U) so that
//              W[i] = Y[i] \cap U.
//           Thus H_P(U) offers the same connectivity information as G_P(U).  
//           ...
///          U consists of those nodes in Y \ X_P that have no parent in Y.

// Question: Should we implement the graph using directed edges?
// Answer:   Not sure.  The connectivity problem does not use directed edges.
//           But it would be easier to check if a node was a root without looking back at the tree
//           if the graph itself has directed edges.  Possibly we could use the 'mark' info instead.


// Question: How do we find the single label for a connected component of size 1?

// Reference: Fast Compatibility Testing for Rooted Phylogenetic Trees
//            Yun Deng, David Fernandex-Baca, Algorithmica (2018) 80:2543-2477
//
//            The overall layout comes from Algorithm 3: BuildST_n(P), lines L1-L23.
//
//            However, in order to implement the steps, we need to
//            initialize (5.2) and maintain (5.3) the data structures mentioned in
//            section (5.1).
//

// Construct a tree that is compatible with all the trees in the profile, and return a null pointer if this is not possible.
//
// Instead of an entire tree, we pass in an internal node from each tree, which is interpreted as a tree where all internal nodes
//   outside of the subtree for that node have been collapsed.  This allows us to avoid creating a new copy of the tree with
//   collapsed nodes.

unique_ptr<Tree_t> BUILD_ST(const vector<const node_t*>& profile)
{
    std::queue<tuple<unique_ptr<connected_component_t>,node_t*>> Q;

    node_t* root = new node_t(nullptr);

//FIXME: info_for_vertex should be replaced with O(1)-lookup vertex attributes.

// L1. Construct display graph H_P(U_init)
    auto [H, Y_inits] = display_graph_from_profile(profile);

// L2.  ENQUEUE(Q, (U_init, null) )
    for(auto& Y_init: Y_inits)
        Q.push({std::move(Y_init), nullptr});

// L3.  while Q is not empty do
    while(not Q.empty())
    {
// L4.  | (U, pred) = DEQUEUE(Q)
        auto [Y,pred] = std::move(Q.front()); Q.pop();

//        assert(Y->count == Y->labels_for_component().size());

// L5.  | Create a node r_U and set parent(r_U) = pred
        node_t* r_U = nullptr;
        if (pred)
        {
            // Create a new node pointing to a higher-level problem.
            r_U = new node_t(pred);
        }
        else
        {
            // This is one of the highest-level problems.
            // Share the root among all the different problems.
            r_U = root;
        }


// L6.  | if |L(U)| = 1 then
        if (Y->count == 1)
        {
// L7.  | | let l be the label in L(U)
            auto labels = Y->labels_for_component();
            assert(labels.size() == 1);
// L8.  | | label r_U with l
            r_U->set_ott_id(labels[0]);
// L9.  | | continue
            continue;
        }

// L10. | if |L(U)| = 2 then
        if (Y->count == 2)
        {
// L11. | | Let l_1, l_2 be the two labels
            auto labels = Y->labels_for_component();
            assert(labels.size() == 2);
// L12. | | foreach j \in [2] do            
            for(auto label: labels)
            {
// L13. | | | Create a node r_j, labelled l_j
// L14. | | | Set parent(r_j) = r_J                
                auto r = new node_t(r_U);
                r->set_ott_id(label);
            }
// L15. | | continue
            continue;
        }

        // The list of semi-universal nodes in U (which corresponds to Y)
        auto SU = Y->semi_universal_nodes_for_position();

        // See Lemma 7.
        //  After replacing each semi-universal node in U with its children,
        //  (i.e. computing the successor of ), there are no semi-universal
        //  nodes in left in U.
        //  BDR: This should get cleared as we modify Y.map(i)
        //    Y->semiU = {};

        // BDR: If we only over split components, then we never need to find
        //      old components to remove them from the list.
        vector<unique_ptr<connected_component_t>> Ws;
        Ws.push_back(std::move(Y));

//      /* Compute the successor of U. */
// L16. | foreach semi-universal node v \in U do
        for(auto parent: SU)
        {
            // Look at recorded component for this vertex.
            auto Y1 = H->vertex_info(parent).get_component();

            int i = *H->vertex_info(parent).tree_index;

            // This is the position U restricted to tree i = U \cap L(i).
            auto& U_i = Y1->marked_vertices_for_tree.at(i);
            assert(U_i.size() == 1);

// L17. | U = (U \ {v}) \cup Ch(v)
            assert(H->in_degree(parent) == 0);
            for(auto [e,e_end]=  H->out_edges(parent); e != e_end; e++)
            {
                auto child = H->target( *e );
                Y1->insert_marked_vertex(child);
            }
            Y1->erase_marked_vertex(parent);

            // By Lemma 7, each semi-universal node always has more than 1
            // child, and so we don't create any new semi-universal nodes.
            // Therefore we do not need to update Y1->semiU.
            auto U_i_tmp = U_i;
            assert(U_i.size() > 1);

            // Remove the edges (parent, child), creating new components, and updating mark, map, and semiU
            for(auto child: U_i_tmp)
                if (H->remove_edge(parent,child))
                {
                    assert(H->vertex_info(child).is_marked());
                    auto Y1 = H->vertex_info(child).get_component();
                    auto Y2 = split_component(Y1, parent, child);
                    if (Y2)
                        Ws.push_back(std::move(Y2));
                }
        }
// L18. Let W_1, W_2, ... W_p be the connected components of H_P(U)

// L19. | if p = 1 then
        assert(Ws.size() > 0);
        if (Ws.size() == 1)
// L20. | | return incompatible
            return {};

// L21. | foreach j \in |p| do
        for(auto& W: Ws)
        {            
// L22. | | ENQUEUE(Q,(W_j, r_U))
            Q.push({std::move(W),r_U});
        }
    }

// L23. return the tree with root r_{U_unit}
    return unique_ptr<Tree_t>{new Tree_t(root)};
}

static vector<int> indices;

/// Construct a tree with all the splits mentioned, and return a null pointer if this is not possible
unique_ptr<Tree_t> BUILD(const vector<int>& tips, const vector<const RSplit*>& splits) {
#pragma clang diagnostic ignored  "-Wsign-conversion"
#pragma clang diagnostic ignored  "-Wsign-compare"
#pragma clang diagnostic ignored  "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored  "-Wsign-compare"
    std::unique_ptr<Tree_t> tree(new Tree_t());
    tree->create_root();
    // 1. First handle trees of size 1 and 2
    if (tips.size() == 1) {
        tree->get_root()->set_ott_id(*tips.begin());
        return tree;
    } else if (tips.size() == 2) {
        auto Node1a = tree->create_child(tree->get_root());
        auto Node1b = tree->create_child(tree->get_root());
        auto it = tips.begin();
        Node1a->set_ott_id(*it++);
        Node1b->set_ott_id(*it++);
        return tree;
    }
    // 2. Initialize the mapping from elements to components
    vector<int> component;       // element index  -> component
    vector<list<int> > elements;  // component -> element indices
    for(int k=0;k<indices.size();k++)
        assert(indices[k] == -1);
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
        for(int id: tips)
            indices[id] = -1;
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
            // indices[i] != -1 checks if x is in the current tip set.
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
        add_subtree(tree->get_root(), *subtree);
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
    for(auto nd: iter_leaf_const(taxonomy)) {
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
    for(auto nd: iter_post_const(summary)) {
        if (nd->get_parent() and nd->get_name().size() and not nd->has_ott_id()) {
            LOG(WARNING)<<"Named taxonomy node has no OTTID!  Not checking for incertae sedis placement.";
            return {};
        }
    }
    map<const Tree_t::node_type*,const Tree_t::node_type*> placements;
    auto node_from_id = get_ottid_to_const_node_map(taxonomy);
    for(auto nd: iter_post_const(summary)) {
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
    for(auto nd: iter_post_const(tree)) {
        // 1. Initialize set for this node.
        auto & inc = include[nd];
        // 2. Add OttId for tip nodes
        if (nd->is_tip()) {
            inc.insert(nd->get_ott_id());
        } else if (nd == tree.get_root()) {
            continue;
        }
        // 3. Add Ids of children only if they are NOT incertae sedis
        for(auto nd2: iter_child_const(*nd)) {
            if (not incertae_sedis.count(nd2->get_ott_id())) {
                auto& inc_child = nd2->get_data().des_ids;
                inc.insert(begin(inc_child),end(inc_child));
            }
        }
    }
    return include;
}

template<typename Tree_T>
map<typename Tree_t::node_type const*, set<OttId>> construct_exclude_sets(const Tree_t& tree, const set<OttId>& incertae_sedis) {
    map<typename Tree_t::node_type const*, set<OttId>> exclude;
    // 1. Set exclude set for root node to the empty set.
    exclude[tree.get_root()];
    for(auto nd: iter_pre_const(tree)) {
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

// % otc-solve-subproblem ott596112.tre
// % otc-solve-subproblem ott497126.tre   OK, used to FAIL is_consistent == is_consistent2
// % otc-solve-subproblem ott5553749.tre
// % otc-solve-subproblem ott56610.tre
// % otc-solve-subproblem ott1041547.tre
// % otc-solve-subproblem ott5551466.tre  OK, used to FAIL is_consistent == is_consistent2

/*
Find test-cases that are not too large:
for i in $(<files-by-size.txt) ; do echo $i ; if ! otc-solve-subproblem $i ; then echo $i >> bad ; fi ; done
 - OK, the problem comes from going through the tree a node at a time.
 - if have ((a1,a2)a,(b1,b2)b), the we test {stuff,a} and keep a, and then {stuff,b} BUT WITHOUT 'a'.
ott791895.tre
ott216628.tre
...
*/

/*
 ott50860.tre

 OK, so we start out with a component with root nodes 0, 17, 34, 48.
 They are in components 1, 18, 18, and 49. (So we start out with three different components.)
 - 0: We replace 0 with its children {1,6} in the frontier.
 */

/*
 * This one is slow: ott1020133.tre
 * This one seems not to terminate: ott1041547.tre
 * - debugging it by gdb/^C/back-trace seems like it could be effective.
 */

template<typename T>
vector<T> append(const vector<T>& v1, const vector<T>& v2)
{
    auto v3 = v1;
    v3.insert(v3.end(), v2.begin(), v2.end());
    return v3;
}

int non_leaf_children(const node_t* start)
{
    int count = 0;
    for(auto node = start->get_first_child(); node ; node = node->get_next_sib())
        if (not node->is_tip())
            count++;
    return count;
}

vector<const node_t*> add_node_remove_children(const vector<const node_t*>& nodes, const node_t* focal_node)
{
    vector<const node_t*> nodes2;
    nodes2.reserve(nodes.size());
    for(auto node: nodes)
        if (node->get_parent() != focal_node)
            nodes2.push_back(node);

    nodes2.push_back(focal_node);
    return nodes2;
}

node_t* copy_subtree(const node_t* subtree)
{
    // 1. Copy all the nodes and record the correspondence.
    std::unordered_map<const node_t*, node_t*> node_map;
    for(auto nd: iter_post_n_const(*subtree))
    {
        auto nd2 = new node_t(nullptr);

        if (nd->has_ott_id())
            nd2->set_ott_id(nd->get_ott_id());

        if (nd->get_name().size())
            nd2->set_name(nd->get_name());

        node_map.insert({nd, nd2});
    }

    // 2. Link corresponding nodes to their corresponding parents
    for(auto& [nd,nd2]: node_map)
    {
        auto p = nd->get_parent();
        auto p2_it = node_map.find(p);

        if (p2_it == node_map.end())
        {
            assert(nd == subtree);
        }
        else
        {
            auto p2 = p2_it->second;
            p2->add_child(nd2);
        }
    }

    // 3. Return the node corresponding to the top of the subtree
    return node_map.at(subtree);
}

unique_ptr<Tree_t> merge_subtrees(const vector<const node_t*>& subtree_nodes)
{
    assert(not subtree_nodes.empty());

    auto root = new node_t(nullptr);

    for(auto node: subtree_nodes)
        root->add_child(copy_subtree(node));

    return unique_ptr<Tree_t>(new Tree_t(root));
}

/// Get the list of splits, and add them one at a time if they are consistent with previous splits
unique_ptr<Tree_t> combine(const vector<unique_ptr<Tree_t>>& trees, const set<OttId>& incertae_sedis, bool verbose) {
    // 0. Standardize names to 0..n-1 for this subproblem
    const auto& taxonomy = trees.back();
    auto all_leaves = taxonomy->get_root()->get_data().des_ids;
    // index -> id
    vector<OttId> ids;
    // id -> index
    map<OttId,int> id_map;
    for(OttId id: all_leaves) {
        int i = ids.size();
        id_map[id] = i;
        ids.push_back(id);
        assert(id_map[ids[i]] == i);
        assert(ids[id_map[id]] == id);
    }
    auto remap = [&id_map](const set<OttId>& argIds) {return remap_ids(argIds, id_map);};
    vector<int> all_leaves_indices;
    for(int i=0;i<all_leaves.size();i++) {
        all_leaves_indices.push_back(i);
    }
    indices.resize(all_leaves.size());
    for(auto& i: indices) {
        i=-1;
    }
    /// Incrementally add splits from @splits_to_try to @consistent if they are consistent with it.
    vector<RSplit> consistent;
    auto add_split_if_consistent = [&](auto nd, RSplit&& split) {
            consistent.push_back(std::move(split));

            auto result = BUILD(all_leaves_indices, consistent);
            if (not result) {
                consistent.pop_back();
                if (verbose and nd->has_ott_id()) {
                    LOG(INFO) << "Reject: ott" << nd->get_ott_id() << "\n";
                }
                return false;
            } else if (verbose and nd->has_ott_id()) {
                LOG(INFO) << "Keep: ott" << nd->get_ott_id() << "\n";
            }
            return true;
        };
    vector<const node_t*> consistent_trees;
    auto remove_split_if_inconsistent = [&](auto focal_node, vector<const node_t*>& consistent_nodes, int nl)
                                        {
                                            // 1. Append the focal node and remove its descendants.
                                            //    There should only be DIRECT descendants on the `consistent_nodes` list.
                                            auto consistent_nodes2 = add_node_remove_children(consistent_nodes, focal_node);

                                            // 2. This is a single tree containing all checked and purified subtrees of the current tree, PLUS the focal node.
                                            auto current_tree = merge_subtrees(consistent_nodes2);
                                            assert(n_leaves(*current_tree) == nl);
                                            // This tree will be destroyed at the end of this routine.

                                            // 3. Check compatibility
                                            consistent_trees.push_back(current_tree->get_root());
                                            auto tree = BUILD_ST( consistent_trees );
                                            consistent_trees.pop_back();

                                            if (not tree)
                                                collapse_split_and_del_node(focal_node);
                                            else
                                                std::swap(consistent_nodes, consistent_nodes2);
                                            return tree;
                                        };
    // 1. Find splits in order of input trees
    vector<Tree_t::node_type const*> compatible_taxa;
    for(int i=0;i<trees.size();i++)
    {
        const auto& tree = trees[i];
        auto root = tree->get_root();
        const auto leafTaxa = root->get_data().des_ids;
        const auto leafTaxaIndices = remap(leafTaxa);

#ifndef NDEBUG
#pragma clang diagnostic ignored  "-Wunreachable-code-loop-increment"
        for(const auto& leaf: set_difference_as_set(leafTaxa, all_leaves)) {
            throw OTCError() << "OTT Id " << leaf << " not in taxonomy!";
        }
#endif
        // Handle the taxonomy tree specially when it has Incertae sedis taxa.
        if (i == trees.size()-1 and not incertae_sedis.empty()) {
            auto exclude = construct_exclude_sets<Tree_t>(*tree, incertae_sedis);

            for(auto nd: iter_post_const(*tree)) {
                if (not nd->is_tip() and nd != root) {
                    // construct split
                    const auto descendants = remap(nd->get_data().des_ids);
                    const auto nondescendants = remap(exclude[nd]);
                    if (add_split_if_consistent(nd, split_from_include_exclude(descendants, nondescendants))) {
                        compatible_taxa.push_back(nd);
                    }
                }
            }
        } else if (i == trees.size()-1) {
            for(auto nd: iter_post_const(*tree)) {
                if (not nd->is_tip() and nd != root) {
                    const auto descendants = remap(nd->get_data().des_ids);
                    if (add_split_if_consistent(nd, RSplit{descendants, leafTaxaIndices})) {
                        compatible_taxa.push_back(nd);
                    }
                }
            }
        } else {
            int consistent_count = 0;

            vector<node_t*> nodes_to_check;
            vector<const node_t*> consistent_nodes;
            // The POST order here means that we always visit a node after all of its children and before any of its ancestors.
            // Therefore we consistent_nodes will always contain all the tips of this tree, if we use add_node_remove_children( ).
            for(auto nd: iter_post(*tree))
            {
                if (nd == root) continue;

                if (nd->is_tip())
                    consistent_nodes.push_back(nd);
                else
                    nodes_to_check.push_back(nd);
            }

            // There are sibling subtrees that have been checked before us.
            // When checking consistency we need to include these.
            int nl = n_leaves(*tree);
            for(auto nd: nodes_to_check)
            {
                const auto descendants = remap(nd->get_data().des_ids);
                bool is_consistent = add_split_if_consistent(nd, RSplit{descendants, leafTaxaIndices});
                LOG(WARNING)<<"Trying partition "<<nd<<" in tree "<<i<<": "<<is_consistent;
                for(auto d: descendants)
                    LOG(DEBUG)<<"    "<<ids[d]<<"  ("<<d<<")";

                bool is_consistent2 = (bool)remove_split_if_inconsistent(nd, consistent_nodes,nl);
                assert(is_consistent == is_consistent2);
                if (is_consistent != is_consistent2)
                    throw OTCError() << "is_consistent != is_consistent2";
                if (is_consistent) consistent_count++;
            }
            if (consistent_count > 0)
                consistent_trees.push_back(tree->get_root());
        }
    }
    // 2. Construct final tree and add names
    auto tree = BUILD(all_leaves_indices, consistent);
    for(auto nd: iter_pre(*tree)) {
        if (nd->is_tip()) {
            int index = nd->get_ott_id();
            nd->set_ott_id(ids[index]);
        }
    }
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
            for(auto nd: iter_pre_const(*tree)) {
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
            for(auto nd: iter_pre_const(*tree)) {
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
        // 7. Perform the synthesis
        auto & taxonomy = *trees.back();
        compute_depth(taxonomy);
        auto tree = combine(trees, incertae_sedis, verbose);
        // 8. Set the root name (if asked)
        // FIXME: This could be avoided if the taxonomy tree in the subproblem always had a name for the root node.
        if (setRootName) {
            tree->get_root()->set_name(args["root-name"].as<string>());
        }
        // 9. Write out the summary tree.
        write_tree_as_newick(std::cout, *tree);
        std::cout << "\n";
        // 10. Find placements
        auto placements = check_placement(*tree, taxonomy);
        for(auto & p: placements) {
            auto placed = p.first;
            auto parent = p.second;
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
