#ifndef OTCETERA_NODE_NAMING_H
#define OTCETERA_NODE_NAMING_H

#include <string>
#include <vector>
#include <unordered_set>
#include "otc/otc_base_includes.h"
#include "otc/tree_operations.h"

namespace otc {

std::string make_name(const std::string& prefix, OttId number);
std::string make_mrca_name(OttId number1, OttId number2);
template<typename N>
OttId smallest_child(const N* node);
template<typename N>
OttId& smallest_child(N* node);
template<typename T>
void calculate_smallest_child(T& tree);
template<typename T>
void sort_by_smallest_child(T& tree);
template<typename T>
void name_unnamed_nodes(const T & tree);

template<typename N>
inline OttId smallest_child(const N * node) {
    return node->get_data().smallest_child;
}

template<typename N>
inline OttId& smallest_child(N * node) {
    return node->get_data().smallest_child;
}

template<typename T>
void calculate_smallest_child(T& tree) {
    for (auto nd: iter_post(tree)) {
        if (nd->is_tip()) {
            smallest_child(nd) = nd->get_ott_id();
        } else {
            auto sc = smallest_child(nd->get_first_child());
            for(auto c: iter_child(*nd)) {
                sc = std::min(sc, smallest_child(c));
            }
            smallest_child(nd) = sc;
        }
    }
}

template<typename T>
void sort_by_smallest_child(T& tree) {
    const std::vector<typename T::node_type*> nodes = all_nodes(tree);
    for (auto nd: nodes) {
        std::vector<typename T::node_type*> children;
        while (nd->has_children()) {
            auto x = nd->get_first_child();
            x->detach_this_node();
            children.push_back(x);
        }
        std::sort(begin(children),
                  end(children),
                  [](const auto& nd1, const auto& nd2)
                  {return smallest_child(nd1) < smallest_child(nd2);});
        while (not children.empty()) {
            auto x = children.back();
            children.pop_back();
            nd->add_child_at_front(x);
        }
    }
}

template<typename T>
void name_unnamed_nodes(T & tree) {
    calculate_smallest_child(tree);
    sort_by_smallest_child(tree);
    
    std::unordered_set<std::string> names;
    for(auto nd:iter_pre(tree)) {
        if (nd->get_name().size()) {
            names.insert(nd->get_name());
        }
    }
    // Remove unnamed nodes w/ no OTT Id that are monotypic.
    std::vector<typename T::node_type*> remove;
    for (auto nd:iter_pre(tree)) {
        if (not nd->has_ott_id() and nd->get_name().empty() and nd->is_outdegree_one_node()){
            remove.push_back(nd);
        }
    }
    for (auto nd: remove) {
        //auto parent = nd->get_parent();
        auto child = nd->get_first_child();
        child->detach_this_node();
        nd->add_sib_on_right(child);
        nd->detach_this_node();
        delete nd;
    }
    
    OttId id = 1;
    for(auto nd:iter_pre(tree)) {
        if (nd->has_ott_id()) {
            nd->set_name("ott"+std::to_string(nd->get_ott_id()));
        } else if (not nd->get_name().size()) {
            assert(not nd->is_tip());
            assert(not nd->is_outdegree_one_node());

            auto id1 = smallest_child(nd->get_first_child());
            auto id2 = smallest_child(nd->get_first_child()->get_next_sib());
            auto name = make_mrca_name(id1,id2);
            if (names.count(name)) {
                throw OTCError()<<"Synthesized name '"<<name<<"' already exists in the tree!";
            }
            nd->set_name(name);
            names.insert(name);
        }
        id++;
    }
}

inline std::string make_name(const std::string& pre, OttId number) {
    return pre + std::to_string(number);
}

inline std::string make_mrca_name(OttId number1, OttId number2) {
    const static std::string mrca_prefix = "mrcaott";
    return mrca_prefix + std::to_string(number1) + "ott" + std::to_string(number2);
}

} //namespace otc
#endif
