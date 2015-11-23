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

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

bool verbose = false;

template <typename T>
std::ostream& operator<<(std::ostream& o, const std::set<T>& s)
{
    auto it = s.begin();
    o<<*it++;
    for(; it != s.end(); it++)
        o<<" "<<*it;
    return o;
}

template <typename T>
std::ostream& operator<<(std::ostream& o, const std::list<T>& s)
{
    auto it = s.begin();
    o<<*it++;
    for(; it != s.end(); it++)
        o<<" "<<*it;
    return o;
}

template <typename T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& s)
{
    auto it = s.begin();
    o<<*it++;
    for(; it != s.end(); it++)
        o<<" "<<*it;
    return o;
}

string newick(const Tree_t &t)
{
    std::ostringstream s;
    writeTreeAsNewick(s, t);
    return s.str();
}

bool get_bool(const string& arg, const string& context="")
{
    if (arg == "true" or arg == "yes" or arg == "True" or arg == "Yes")
        return true;
    else if (arg == "false" or arg == "no" or arg == "False" or arg == "No")
        return false;
    else
        throw OTCError()<<context<<"'"<<arg<<"' is not a recognized boolean value.";
}

bool chopRoot = false;
bool handleChopRoot(OTCLI & otCLI, const std::string &)
{
    chopRoot = true;
    return true;
}

string prefix = "node";
bool handlePrefix(OTCLI & otCLI, const std::string & arg)
{
    prefix = arg;
    return true;
}

string makeName(const string& prefix, int number)
{
    return prefix + std::to_string(number);
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-graft-solutions",
                "Takes a series of tree files, which are treated as subproblem solutions.\n"
                "Each solution tree should have an OTT Id at the root.\n",
                "solutions.tre");

    otCLI.addFlag('p',
                  "Prefix for unnamed nodes",
                  handlePrefix,
                  true);
    
    otCLI.addFlag('c',
                  "Chop of the root node",
                  handleChopRoot,
                  false);
    

    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};

    if (argc < 2)
        throw OTCError("No trees provided!");

    // I think multiple subproblem files are essentially concatenated.
    // Is it possible to read a single subproblem from cin?
    if (treeProcessingMain<Tree_t>(otCLI, argc, argv, get, nullptr, 1))
        std::exit(1);

    if (chopRoot)
    {
        Tree_t& tree = *trees[0];
        auto newRoot = tree.getRoot()->getFirstChild();
        newRoot->detachThisNode();
        tree._setRoot(newRoot);
    }

    verbose = otCLI.verbose;

    if (trees.empty())
        throw OTCError("No trees loaded!");

    // Add names to unnamed nodes
    for(const auto& tree: trees)
    {
        std::unordered_set<string> names;
        for(auto nd:iter_pre(*tree))
            if (nd->getName().size())
                names.insert(nd->getName());

        int id = 1;
        for(auto nd:iter_pre(*tree))
        {
            if (nd->hasOttId())
                nd->setName("ott"+std::to_string(nd->getOttId()));
            else if (not nd->getName().size())
            {
                string name = makeName(prefix,id);
                if (names.count(name))
                    throw OTCError()<<"Synthesized name '"<<name<<"' already exists in the tree!";
                nd->setName(makeName(prefix,id));
            }
            id++;
        }
    }

    for(const auto& tree: trees)
    {
        writeTreeAsNewick(std::cout, *tree);
        std::cout<<"\n";
    }
}
