#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <iterator>
#include <sstream>

#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_iter.h"
using namespace otc;
using std::vector;
using std::unique_ptr;
using std::set;
using std::list;
using std::map;
using std::string;
using namespace otc;

struct RTNodeSmallestChild
{
    int smallestChild = 0;
};

using Tree_t = RootedTree<RTNodeSmallestChild, RTreeNoData>;

inline int smallestChild(const Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

inline int& smallestChild(Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

static std::string mrca_prefix = "mrcaott";

string makeName(const string& prefix, int number);

string makeName(const string& pre, int number) {
    return pre + std::to_string(number);
}

string makeMRCAName(int number1, int number2) {
    return mrca_prefix + std::to_string(number1) + "ott" + std::to_string(number2);
}

void calculateSmallestChild(Tree_t& T)
{
    for(auto nd: iter_post(T))
        if (nd->isTip())
            smallestChild(nd) = nd->getOttId();
        else
        {
            int sc = smallestChild(nd->getFirstChild());
            for(auto c: iter_child(*nd))
                sc = std::min(sc, smallestChild(c));
            smallestChild(nd) = sc;
        }
}

void sortBySmallestChild(Tree_t& T)
{
    vector<Tree_t::node_type*> nodes;
    for(auto nd: iter_post(T))
        if (not nd->isTip())
            nodes.push_back(nd);

    for(auto nd: nodes)
    {
        vector<Tree_t::node_type*> children;
        while(nd->hasChildren())
        {
            auto x = nd->getFirstChild();
            x->detachThisNode();
            children.push_back(x);
        }
        std::sort( begin(children), end(children),
                   [](const auto& nd1, const auto& nd2)
                   {return smallestChild(nd1) < smallestChild(nd2);}
            );
        while(not children.empty())
        {
            auto x = children.back();
            children.pop_back();
            nd->addChildAtFront(x);
        }
    }
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-name-unnamed-nodes",
                "Takes a series of tree files and writes each tree with names computed for unnamed nodes.\n"
                "Nodes with OTT ids are written as ott#####, with longer names suppressed.\n",
                "tree1.tre [tree2.tre ... treeN.tree]");
    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};
    if (argc < 2) {
        throw OTCError("No trees provided!");
    }
    if (treeProcessingMain<Tree_t>(otCLI, argc, argv, get, nullptr, 1)) {
        return 1;
    }
    if (trees.empty()) {
        throw OTCError("No trees loaded!");
    }
    // Add names to unnamed nodes
    for(const auto& tree: trees) {
        calculateSmallestChild(*tree);
        sortBySmallestChild(*tree);
        
        std::unordered_set<string> names;
        for(auto nd:iter_pre(*tree)) {
            if (nd->getName().size()) {
                names.insert(nd->getName());
            }
        }

        // Remove unnamed nodes w/ no OTT Id that are monotypic.
        vector<Tree_t::node_type*> remove;
        for(auto nd:iter_pre(*tree))
            if (not nd->hasOttId() and nd->getName().empty() and nd->isOutDegreeOneNode())
                remove.push_back(nd);
        for(auto nd: remove)
        {
            auto parent = nd->getParent();
            auto child = nd->getFirstChild();
            child->detachThisNode();
            nd->addSibOnRight(child);
            nd->detachThisNode();
            delete nd;
        }
        
        int id = 1;
        for(auto nd:iter_pre(*tree)) {
            if (nd->hasOttId()) {
                nd->setName("ott"+std::to_string(nd->getOttId()));
            } else if (not nd->getName().size()) {
                assert(not nd->isTip());
                assert(not nd->isOutDegreeOneNode());

                int id1 = smallestChild(nd->getFirstChild());
                int id2 = smallestChild(nd->getFirstChild()->getNextSib());
                string name = makeMRCAName(id1,id2);
                if (names.count(name)) {
                    throw OTCError()<<"Synthesized name '"<<name<<"' already exists in the tree!";
                }
                nd->setName(name);
                names.insert(name);
            }
            id++;
        }
    }
    for(const auto& tree: trees) {
        writeTreeAsNewick(std::cout, *tree);
        std::cout<<"\n";
    }
    return 0;
}
