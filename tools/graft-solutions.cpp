#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <iterator>

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

/// Create a mapping from name -> id
map<string, long> createIdsFromNames(const vector<unique_ptr<Tree_t>>& trees)
{
    long id = 1;
    map<string,long> name_to_id;
    for(const auto& tree: trees)
        for(auto nd: iter_post_const(*tree))
            if (nd->getName().size())
            {
                string name = nd->getName();
                auto it = name_to_id.find(name);
                if (it == name_to_id.end())
                    name_to_id[name] = id++;
            }
            else if (nd->isTip())
                throw OTCError()<<"Tip has no label!";

    return name_to_id;
}  

/// Set ids on the tree based on the name
void setIdsFromNames(Tree_t& tree, const map<string,long>& name_to_id)
{
    for(auto nd: iter_post(tree))
        if (nd->getName().size())
        {
            string name = nd->getName();
            auto it = name_to_id.find(name);
            if (it == name_to_id.end())
                throw OTCError()<<"Can't find label '"<<name<<"' in taxonomy!";
            auto id = it->second;
            nd->setOttId(id);
        }
        else if (nd->isTip())
            throw OTCError()<<"Tree tip has no label!";
}

string addOttId(const string s, long id)
{
    string tag = "ott" + std::to_string(id);
    if (not s.size())
        return tag;
    else
        return s + " " + tag;
}

void relabelWithOttId(Tree_t& T)
{
    for(auto nd: iter_pre(T))
        if (nd->hasOttId())
            nd->setName(addOttId(nd->getName(),nd->getOttId()));
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

bool handleRequireOttIds(OTCLI & otCLI, const std::string & arg)
{
    otCLI.getParsingRules().setOttIds = get_bool(arg,"-o: ");
    return true;
}

string rootName = "";

bool handleRootName(OTCLI& otCLI, const std::string & arg)
{
    rootName = arg;
    return true;
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-graft-solutions",
                "Takes a series of tree files, which are treated as subproblem solutions.\n"
                "Each solution tree should have an OTT Id at the root.\n",
                "solutions.tre");

    otCLI.addFlag('o',
                  "Require OTT ids.  Defaults to true",
                  handleRequireOttIds,
                  true);

    otCLI.addFlag('n',
                  "Rename the root to this name",
                  handleRootName,
                  true);

    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};

    if (argc < 2)
        throw OTCError("No solutions provided!");

    // I think multiple subproblem files are essentially concatenated.
    // Is it possible to read a single subproblem from cin?
    if (treeProcessingMain<Tree_t>(otCLI, argc, argv, get, nullptr, 1))
        std::exit(1);

    verbose = otCLI.verbose;

    if (trees.empty())
        throw OTCError("No trees loaded!");

    bool setOttIds = otCLI.getParsingRules().setOttIds;

    vector<int> no_root_label;
    for(int i=0;i<trees.size();i++)
        if (trees[i]->getRoot()->getName().empty())
            no_root_label.push_back(i);

    if (no_root_label.size() > 1)
    {
        OTCError e;
        e<<no_root_label.size()<<" trees have an unlabelled root!\n";
        int n = std::min(10,int(no_root_label.size()));
        e<<"  They are trees "<<no_root_label[0];
        for(int i=1;i<n;i++)
            e<<", "<<no_root_label[i];
        if (n > 10)
            e<<" ...";
        else
            e<<".";
        throw e;
    }

    if (not setOttIds)
    {
        auto name_to_id = createIdsFromNames(trees);
        for(auto& tree: trees)
            setIdsFromNames(*tree, name_to_id);
    }

    std::unordered_map<long,Tree_t::node_type*> my_leaf;

    // Find the nodes where we would like to graft a tree.
    // Each tip id should occur only once as a tip.
    for(const auto& tree: trees)
        for(auto nd:iter_pre(*tree))
            if (nd->isTip())
            {
                assert(nd->hasOttId());
                long id = nd->getOttId();
                if (my_leaf.find(id) != my_leaf.end())
                {
                    if (setOttIds)
                        throw OTCError()<<"OTT Id "<<id<<" occurs at multiple tips!";
                    else
                        throw OTCError()<<"Label '"<<nd->getName()<<"' occurs at multiple tips!";
                }
                my_leaf[id] = nd;
            }

    // Check that we don't have multiple examples of the same subproblem.
    // Each root id should occur only once as a root.
    std::unordered_set<long> root_ids;
    for(int i=0;i<trees.size();i++)
    {
        auto root = trees[i]->getRoot();
        long id = root->getOttId();
        if (not root_ids.count(id))
            root_ids.insert(id);
        else
            if (setOttIds)
                throw OTCError()<<"OTT Id "<<id<<" occurs at the root of multiple trees!";
            else
                throw OTCError()<<"Label '"<<root->getName()<<"' occurs at the root of multiple trees!";
    }
    
    // Glue each root into its corresponding tip.
    vector<unique_ptr<Tree_t>> roots;
    for(int i=0;i<trees.size();i++)
    {
        long id = trees[i]->getRoot()->getOttId();
        if (not my_leaf.count(id))
        {
            if (otCLI.verbose)
                LOG(INFO)<<"OTT Id "<<id<<" is not a leaf in any subproblem.  Must be a root.\n";
            roots.push_back({});
            std::swap(roots.back(), trees[i]);
        }
        else
        {
            auto nd = my_leaf[id];
            replaceWithSubtree<Tree_t>(nd, *trees[i]);
        }
    }

    if (roots.size() == 1 and not rootName.empty())
        roots[0]->getRoot()->setName(rootName);
    
    for(const auto& tree: roots)
    {
        writeTreeAsNewick(std::cout, *tree);
        std::cout<<"\n";
    }

    if (roots.size() != 1)
        return 1;
}
