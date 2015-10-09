#include <algorithm>
#include <set>
#include <unordered_map>
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

typedef TreeMappedWithSplits Tree_t;

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
      tree.getData().ottIdToNode[id] = nd;
    }
    else if (nd->isTip())
      throw OTCError()<<"Tree tip has no label!";
  
  clearAndfillDesIdSets(tree);
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
  otCLI.getParsingRules().requireOttIds = get_bool(arg,"-o: ");
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
                "Each subproblem tree should have an OTT Id at the root.\n",
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

    bool requireOttIds = otCLI.getParsingRules().requireOttIds;

    for(int i=0;i<trees.size();i++)
      if (trees[i]->getRoot()->getName().empty())
	throw OTCError()<<"Tree "<<i+1<<" has an unlabelled root!";
      else if (requireOttIds and not trees[i]->getRoot()->hasOttId())
	throw OTCError()<<"Tree "<<i+1<<" has no OTT Id for the root!";

    if (not requireOttIds)
    {
      auto name_to_id = createIdsFromNames(trees);
      for(auto& tree: trees)
	setIdsFromNames(*tree, name_to_id);
    }

    std::unordered_map<long,Tree_t::node_type*> my_leaf;

    // Find the nodes where we would like to graft a tree
    for(const auto& tree: trees)
      for(auto nd:iter_pre(*tree))
	if (nd->isTip())
	{
	  assert(nd->hasOttId());
	  long id = nd->getOttId();
	  assert(my_leaf.find(id) == my_leaf.end());
	  my_leaf[id] = nd;
	}


    vector<unique_ptr<Tree_t>> roots;

    for(int i=0;i<trees.size();i++)
    {
      long id = trees[i]->getRoot()->getOttId();
      if (not my_leaf.count(id))
      {
	if (otCLI.verbose)
	  LOG(INFO)<<"OTT Id "<<id<<" is not a leaf in any subproblem.  Must be a root.\n";
	roots.push_back({});
	std::swap(roots.back(),trees[i]);
      }
      else
      {
	auto nd = my_leaf[id];
	replaceWithSubtree<Tree_t>(nd, *trees[i]);
      }
    }
    for(const auto& tree: roots)
    {
      writeTreeAsNewick(std::cout, *tree);
      std::cout<<"\n";
    }
}
