#ifndef OTCETERA_NODE_NAMING_H
#define OTCETERA_NODE_NAMING_H

#include <string>
#include <vector>
#include <unordered_set>
#include "otc/otc_base_includes.h"
#include "otc/tree_operations.h"

namespace otc {

std::string makeName(const std::string& prefix, long number);
std::string makeMRCAName(long number1, long number2);
template<typename N>
OttId smallest_child(const N* node);
template<typename N>
OttId& smallest_child(N* node);
template<typename T>
void calculateSmallestChild(T& tree);
template<typename T>
void sortBySmallestChild(T& tree);
template<typename T>
void nameUnamedNodes(const T & tree);

template<typename N>
inline OttId smallest_child(const N * node) {
    return node->get_data().smallest_child;
}

template<typename N>
inline OttId& smallest_child(N * node) {
    return node->get_data().smallest_child;
}

template<typename T>
void calculateSmallestChild(T& tree) {
    for (auto nd: iter_post(tree)) {
        if (nd->isTip()) {
            smallest_child(nd) = nd->get_ott_id();
        } else {
            auto sc = smallest_child(nd->getFirstChild());
            for(auto c: iter_child(*nd)) {
                sc = std::min(sc, smallest_child(c));
            }
            smallest_child(nd) = sc;
        }
    }
}

template<typename T>
void sortBySmallestChild(T& tree) {
    const std::vector<typename T::node_type*> nodes = all_nodes(tree);
    for (auto nd: nodes) {
        std::vector<typename T::node_type*> children;
        while (nd->hasChildren()) {
            auto x = nd->getFirstChild();
            x->detachThisNode();
            children.push_back(x);
        }
        std::sort(begin(children),
                  end(children),
                  [](const auto& nd1, const auto& nd2)
                  {return smallest_child(nd1) < smallest_child(nd2);});
        while (not children.empty()) {
            auto x = children.back();
            children.pop_back();
            nd->addChildAtFront(x);
        }
    }
}

template<typename T>
void nameUnamedNodes(T & tree) {
    calculateSmallestChild(tree);
    sortBySmallestChild(tree);
    
    std::unordered_set<std::string> names;
    for(auto nd:iter_pre(tree)) {
        if (nd->get_name().size()) {
            names.insert(nd->get_name());
        }
    }
    // Remove unnamed nodes w/ no OTT Id that are monotypic.
    std::vector<typename T::node_type*> remove;
    for (auto nd:iter_pre(tree)) {
        if (not nd->has_ott_id() and nd->get_name().empty() and nd->isOutDegreeOneNode()){
            remove.push_back(nd);
        }
    }
    for (auto nd: remove) {
        //auto parent = nd->getParent();
        auto child = nd->getFirstChild();
        child->detachThisNode();
        nd->addSibOnRight(child);
        nd->detachThisNode();
        delete nd;
    }
    
    long id = 1;
    for(auto nd:iter_pre(tree)) {
        if (nd->has_ott_id()) {
            nd->setName("ott"+std::to_string(nd->get_ott_id()));
        } else if (not nd->get_name().size()) {
            assert(not nd->isTip());
            assert(not nd->isOutDegreeOneNode());

            auto id1 = smallest_child(nd->getFirstChild());
            auto id2 = smallest_child(nd->getFirstChild()->getNextSib());
            auto name = makeMRCAName(id1,id2);
            if (names.count(name)) {
                throw OTCError()<<"Synthesized name '"<<name<<"' already exists in the tree!";
            }
            nd->setName(name);
            names.insert(name);
        }
        id++;
    }
}

inline std::string makeName(const std::string& pre, long number) {
    return pre + std::to_string(number);
}

inline std::string makeMRCAName(long number1, long number2) {
    const static std::string mrca_prefix = "mrcaott";
    return mrca_prefix + std::to_string(number1) + "ott" + std::to_string(number2);
}

} //namespace otc
#endif
