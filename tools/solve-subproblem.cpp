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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <optional>

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

class dynamic_graph
{
    Graph G;

    map<int,list<Vertex>> vertices_for_component_;

    // This should really be encoded ON the vertex
    map<Vertex,int> component_for_vertex_;

public:
    int n_components() const {return vertices_for_component_.size();}

    int component_for_vertex(Vertex v) const {return component_for_vertex_.at(v);}

    const list<Vertex>& vertices_for_component(int c) const {return vertices_for_component_.at(c);}

          list<Vertex>& vertices_for_component(int c)       {return vertices_for_component_.at(c);}

    int size_of_component(int c) const {return vertices_for_component(c).size();}

    auto in_edges(Vertex v) const {return boost::in_edges(v,G);}

    auto out_edges(Vertex v) const {return boost::out_edges(v,G);}

    Vertex source(Edge e) const {return boost::source(e,G);}

    Vertex target(Edge e) const {return boost::target(e,G);}

    int out_degree(Vertex v) const
    {
        int count = 0;
        for(auto [e, e_end] = out_edges(v); e != e_end; e++)
            count++;
        return count;
    }

    bool is_internal_node(Vertex v) const
    {
        auto [e, e_end] = out_edges(v);
        return (e != e_end);
    }

    Vertex add_vertex()
    {
        auto v = boost::add_vertex(G);
        int c = n_components();
        vertices_for_component_.insert({c,{v}});
        return v;
    }

    auto add_edge(Vertex u, Vertex v)
    {
        auto e = boost::add_edge(u,v,G);

        int cu = component_for_vertex(u);
        int cv = component_for_vertex(v);

        if (cu != cv)
        {
            if (size_of_component(cu) > size_of_component(cv)) std::swap(cu,cv);

            // Now cu is smaller, or equal.
            auto& lu = vertices_for_component(cu);
            auto& lv = vertices_for_component(cv);
            for(auto uu: lu)
                component_for_vertex_[uu] = cv;
            lv.splice(lv.end(), lu);
            vertices_for_component_.erase(cu);
        }
        return e;
    }

    bool remove_edge(Vertex u, Vertex v)
    {
        int c1 = component_for_vertex(u);
        boost::remove_edge(u,v,G);
        vector<Vertex> from_u;
        vector<Vertex> from_v;
        map<Vertex,int> seen;

        auto try_add_u = [&](Vertex uu)
                             {
                                 auto u_it = seen.find(uu);
                                 if (u_it == seen.end())
                                 {
                                     seen[uu] = 0;
                                     from_u.push_back(uu);
                                     return false;
                                 }
                                 else if (u_it->second == 1)
                                     return true;
                                 else
                                     return false;
                             };

        auto try_add_v = [&](Vertex vv)
                             {
                                 auto v_it = seen.find(vv);
                                 if (v_it == seen.end())
                                 {
                                     seen[vv] = 1;
                                     from_v.push_back(vv);
                                     return false;
                                 }
                                 else if (v_it->second == 0)
                                     return true;
                                 else
                                     return false;
                             };

        try_add_u(u);
        try_add_v(v);

        int i=0,j=0;
        while(i < from_u.size() and j < from_v.size())
        {
            // Check 1 entry from u
            auto uu = from_u[i++];
            for(auto [e, e_end] = in_edges(uu); e != e_end; e++)
                if (try_add_u( source( *e ) )) return false;
            for(auto [e, e_end] = out_edges(uu); e != e_end; e++)
                if (try_add_u( target( *e ) )) return false;

            // Check 1 entry from v
            auto vv = from_v[j++];
            for(auto [e, e_end] = in_edges(vv); e != e_end; e++)
                if (try_add_v( source( *e ) )) return false;
            for(auto [e, e_end] = out_edges(vv); e != e_end; e++)
                if (try_add_v( target( *e ) )) return false;
        }

        if (i >= from_u.size())
            std::swap(from_u, from_v);

        int c2 = n_components();
        for(auto uu: from_u)
            component_for_vertex_[uu] = c2;
        auto& nodes1 = vertices_for_component(c1);
        list<Vertex> nodes2;
        for(auto it = nodes1.begin(); it != nodes1.end();)
        {
            auto next = it; next++;
            if (component_for_vertex_[*it] == c2)
                nodes2.splice(nodes2.end(), nodes1, it, next);
            it = next;
        }
        vertices_for_component_[c2] = std::move(nodes2);

        return true;
    }
};


struct connected_component_t;

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

    bool is_marked() const {return component;}
    void unmark_node() {component = nullptr;}
    void mark_node(connected_component_t* c) {component = c;}

    connected_component_t* get_component() const {assert(component); return component;}
};

bool semi_universal_position(const dynamic_graph& G, const set<Vertex>& vs)
{
    if (vs.size() == 1)
    {
        auto v = *vs.begin();
        if (G.is_internal_node(v)) return true;
    }
    return false;
}

struct connected_component_t
{
    dynamic_graph* G;

    dynamic_graph* get_graph() {return G;}

    // Y.count: the cardinality of Y \cap X_P
    // That is, the number of tip labels in the component.
    int count = -1;

    // Y.map: {(i,L[i])| i is a tree index, and L[i] is not empty}
    //        where L[i] is the set of marked vertices in tree T[i].
    map<int,set<Vertex>> marked_vertices_for_tree;

    void erase_marked_vertex_for_tree(Vertex v, int i)
    {
        auto& nodes = marked_vertices_for_tree.at(i);
        if (semi_universal_position(*G, nodes))
            semiU.erase(i);
        nodes.erase(v);
        if (semi_universal_position(*G, nodes))
            semiU.insert(i);
    }

    void insert_marked_vertex_for_tree(Vertex v, int i)
    {
        auto& nodes = marked_vertices_for_tree.at(i);
        if (semi_universal_position(*G, nodes))
            semiU.erase(i);
        nodes.insert(v);
        if (semi_universal_position(*G, nodes))
            semiU.insert(i);
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

    vector<OttId> labels_for_component(const map<Vertex,vertex_info_t>& info_for_vertex) const
    {
        vector<OttId> labels;
        auto vs = vertices_in_component();
        for(auto v: vs)
            if (auto l = info_for_vertex.at(v).label)
                labels.push_back(*l);
        return labels;
    }

    connected_component_t(dynamic_graph* g):G(g) {}
};

typedef vector<const node_t*> position_t;



tuple<unique_ptr<dynamic_graph>, map<Vertex, vertex_info_t>,unique_ptr<connected_component_t>>
display_graph_from_profile(const vector<const node_t*>& profile)
{
    unique_ptr<dynamic_graph> H;
    map<OttId,Vertex> Labels;
    map<const node_t*,Vertex> node_to_vertex;
    map<Vertex,vertex_info_t> info_for_vertex;
    unique_ptr<connected_component_t> Y_init(new connected_component_t(H.get()));

    for(int i=0; i< profile.size(); i++)
    {
        auto root = profile[i];

        //----------- 5.2 -------------//
        info_for_vertex[node_to_vertex[root]].mark_node(Y_init.get());

        // Walk nodes below root in pre-order (parent before child) so that we can connect children to parents.
        for(auto nd: iter_pre_n_const(root))
        {
            // Add vertex for child node.
            auto v = H->add_vertex();
            node_to_vertex.insert({nd,v});
            info_for_vertex[v].tree_index = i;

            // Add edge from parent node to child node
            if (nd->get_parent())
            {
                auto u = node_to_vertex.at(nd->get_parent());
                H->add_edge(u,v);
            }

            // Add an edge from the node for the label to the tip
            if (nd->is_tip())
            {
                assert(nd->has_ott_id());
                auto id = nd->get_ott_id();

                if (not Labels.count(id))
                {
                    auto label = H->add_vertex();
                    Labels.insert({id,label});
                    info_for_vertex[label].label = id;
                }

                H->add_edge(v,Labels.at(id));
            }
        }
    }

    // * Y_init.count = |L(P)|
    Y_init->count = Labels.size();

    // * Y_init.map consists of all pairs (i, {r(T[i])}), for each i \in [k] ([0..k-1] for our purposes)
    for(int i=0;i<profile.size();i++)
    {
        auto root = profile[i];
        Y_init->marked_vertices_for_tree.insert({i,{node_to_vertex[root]}});
    }

    // Y_init.semiU = [k] ([0..k-1] for our purposes)
    for(int i=0;i<profile.size();i++)
        Y_init->semiU.insert(i);

    return {std::move(H), info_for_vertex, std::move(Y_init)};
}


// Walk the component Y2 containing u and move marked nodes from Y1 to Y2
void split_component(dynamic_graph& H, map<Vertex,vertex_info_t>& info_for_vertex, connected_component_t* Y1, connected_component_t* Y2, Vertex u)
{
    int cu = H.component_for_vertex(u);

    for(auto uu: H.vertices_for_component(cu))
    {
        auto& info = info_for_vertex[uu];
        if (not info.tree_index)
        {
            Y1->count--;
            Y2->count++;
            continue;
        }
        int i = *info.tree_index;
        if (info.is_marked())
        {
            Y2->erase_marked_vertex_for_tree(uu, i);

            Y1->insert_marked_vertex_for_tree(uu, i);
            info.mark_node(Y1);
        }
    }
}

// We need to return a new component corresponding to the second vertex u
unique_ptr<connected_component_t> split_component(dynamic_graph& H, map<Vertex,vertex_info_t>& info_for_vertex, connected_component_t* Y1, Vertex v, Vertex u)
{
    unique_ptr<connected_component_t> W(new connected_component_t(&H));
    auto Y2 = W.get();

    int cv = H.component_for_vertex(v);
    int cu = H.component_for_vertex(u);

    if (H.size_of_component(cu) > H.size_of_component(cv))
        split_component(H, info_for_vertex, Y2, Y1, v);
    else
        split_component(H, info_for_vertex, Y1, Y2, u);

    return W;
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
// Construct a tree that is compatible with all the trees in the profile, and return a null pointer if this is not possible
unique_ptr<Tree_t> BUILD_ST(const vector<const node_t*>& profile)
{
    std::queue<tuple<unique_ptr<connected_component_t>,node_t*>> Q;
    node_t* r_U_init = nullptr;

//FIXME: info_for_vertex should be replaced with O(1)-lookup vertex attributes.

// L1. Construct display graph H_P(U_init)
    auto [H, info_for_vertex, Y_init] = display_graph_from_profile(profile);

// L2.  ENQUEUE(Q, (U_init, null) )
    Q.push({std::move(Y_init), nullptr});

// L3.  while Q is not empty do
    while(not Q.empty())
    {
// L4.  | (U, pred) = DEQUEUE(Q)
        auto [Y,pred] = std::move(Q.front()); Q.pop();

// L5.  | Create a node r_U and set parent(r_U) = pred
        auto r_U = new node_t(pred);

        if (not pred) r_U_init = r_U; // Is there a more elegant way to save r_U_init?

// L6.  | if |L(U)| = 1 then
        if (Y->count == 1)
        {
// L7.  | | let l be the label in L(U)
            auto labels = Y->labels_for_component(info_for_vertex);
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
            auto labels = Y->labels_for_component(info_for_vertex);
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
        Y->semiU = {};

        // BDR: If we only over split components, then we never need to find
        //      old components to remove them from the list.
        vector<unique_ptr<connected_component_t>> Ws;
        Ws.push_back(std::move(Y));

//      /* Compute the successor of U. */
// L16. | foreach semi-universal node v \in U do
        for(auto v: SU)
        {
            // Look at recorded component for this vertex.
            auto Y1 = info_for_vertex[v].get_component();

            int i = *info_for_vertex[v].tree_index;

            // This is the position U restricted to tree i = U \cap L(i).
            auto& U_i = Y1->marked_vertices_for_tree.at(i);
            assert(U_i.size() == 1);

// L17. | U = (U \ {v}) \cup Ch(v)
            U_i.clear();
            info_for_vertex[v].unmark_node();
            for(auto [e,e_end]=  H->out_edges(v); e != e_end; e++)
            {
                auto u = H->target( *e );
                U_i.insert(u);
                info_for_vertex[u].mark_node(Y1);
            }

            // By Lemma 7, each semi-universal node always has more than 1
            // child, and so we don't create any new semi-universal nodes.
            // Therefore we do not need to update Y1->semiU.
            assert(U_i.size() > 1);

            // Remove the edges (v,u), creating new components, and updating mark, map, and semiU
            for(auto u: U_i)
                if (H->remove_edge(v,u))
                    Ws.push_back(split_component(*H, info_for_vertex, Y1, v, u));
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
    return unique_ptr<Tree_t>{new Tree_t(r_U_init)};
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
    auto add_split_if_consistent = [&all_leaves_indices,verbose,&consistent](auto nd, RSplit&& split) {
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
    // 1. Find splits in order of input trees
    vector<Tree_t::node_type const*> compatible_taxa;
    for(int i=0;i<trees.size();i++) {
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
            for(auto nd: iter_post_const(*tree)) {
                if (not nd->is_tip() and nd != root) {
                    const auto descendants = remap(nd->get_data().des_ids);
                    add_split_if_consistent(nd, RSplit{descendants, leafTaxaIndices});
                }
            }
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

int main(int argc, char *argv[]) {
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
