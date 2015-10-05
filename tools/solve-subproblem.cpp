#include <algorithm>
#include <set>
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

struct RSplit
{
  set<int> in;
  set<int> out;
  set<int> all;
  RSplit() = default;
  RSplit(const set<int>& i, const set<int>& a)
    :in(i),all(a)
  {
    out = set_difference_as_set(all,in);
    assert(in.size() + out.size() == all.size());
  }
};

std::ostream& operator<<(std::ostream& o, const RSplit& s)
{
  o<<s.in<<" | "<<s.out;
  return o;
}

/// Merge components c1 and c2 and return the component name that survived
int merge_components(int c1, int c2, vector<int>& component, vector<list<int>>& elements)
{
  if (elements[c2].size() > elements[c1].size())
    std::swap(c1,c2);

  for(int i:elements[c2])
    component[i] = c1;

  elements[c1].splice(elements[c1].end(), elements[c2]);
  
  return c1;
}

bool empty_intersection(const set<int>& x, const set<int>& y)
{
  std::set<int>::const_iterator i = x.begin();
  std::set<int>::const_iterator j = y.begin();
  while (i != x.end() && j != y.end())
  {
    if (*i == *j)
      return false;
    else if (*i < *j)
      ++i;
    else
      ++j;
  }
  return true;
}

bool empty_intersection(const set<int>& xs, const vector<int>& ys)
{
  for(int y: ys)
    if (xs.count(y))
      return false;
  return true;
}

vector<int> indices;

/// Construct a tree with all the splits mentioned, and return a null pointer if this is not possible
unique_ptr<Tree_t> BUILD(const vector<int>& tips, const vector<const RSplit*>& splits)
{
  std::unique_ptr<Tree_t> tree(new Tree_t());
  tree->createRoot();

  // 1. First handle trees of size 1 and 2
  if (tips.size() == 1)
  {
    tree->getRoot()->setOttId(*tips.begin());
    return tree;
  }
  else if (tips.size() == 2)
  {
    auto Node1a = tree->createChild(tree->getRoot());
    auto Node1b = tree->createChild(tree->getRoot());
    auto it = tips.begin();
    Node1a->setOttId(*it++);
    Node1b->setOttId(*it++);
    return tree;
  }

  // 2. Initialize the mapping from elements to components
  vector<int> component;       // element index  -> component
  vector<list<int>> elements;  // component -> element indices
  for(int i=0;i<tips.size();i++)
  {
    indices[tips[i]] = i;
    component.push_back(i);
    elements.push_back({i});
  }

  // 3. For each split, all the leaves in the include group must be in the same component
  for(const auto& split: splits)
  {
    int c1 = -1;
    for(int i: split->in)
    {
      int j = indices[i];
      int c2 = component[j];
      if (c1 != -1 and c1 != c2)
	merge_components(c1,c2,component,elements);
      c1 = component[j];
    }
  }

  // 4. If we can't subdivide the leaves in any way, then the splits are not consistent, so return failure
  if (elements[component[0]].size() == tips.size())
    return {};

  // 5. Create the set of tips in each connected component 
  map<int,vector<int>> subtips;
  for(int c=0;c<tips.size();c++)
  {
    if (c != component[c]) continue;
    
    vector<int>& s = subtips[c];
    for(int i: elements[c])
      s.push_back(tips[i]);
  }

  // 6. Determine the splits that are not satisfied yet and go into each component
  map<int,vector<const RSplit*>> subsplits;
  for(const auto& split: splits)
  {
    int first = indices[*split->in.begin()];
    assert(first >= 0);
    int c = component[first];

    // if none of the exclude group are in the component, then the split is satisfied by the top-level partition.
    if (empty_intersection(split->out, subtips[c])) continue;

    subsplits[c].push_back(split);
  }
  
  // Clear our map from id -> index
  for(int id: tips)
    indices[id] = -1;
  
  // 7. Recursively solve the sub-problems of the partition components
  for(const auto& x: subtips)
  {
    auto subtree = BUILD(x.second, subsplits[x.first]);
    if (not subtree) return {};

    tree->addSubtree(tree->getRoot(), *subtree);
  }

  return tree;
}

/// Copy node names from taxonomy to tree based on ott ids, and copy the root name also
void add_names(unique_ptr<Tree_t>& tree, const unique_ptr<Tree_t>& taxonomy)
{
  clearAndfillDesIdSets(*tree);

  for(auto n1: iter_post(*tree))
    for(auto n2: iter_post(*taxonomy))
      if (n1->getData().desIds == n2->getData().desIds)
	n1->setName( n2->getName());
}

set<int> remap_ids(const set<long>& s1, const map<long,int>& id_map)
{
  set<int> s2;
  for(auto x: s1)
  {
    auto it = id_map.find(x);
    assert(it != id_map.end());
    s2.insert(it->second);
  }
  return s2;
}


/// Get the list of splits, and add them one at a time if they are consistent with previous splits
unique_ptr<Tree_t> combine(const vector<unique_ptr<Tree_t>>& trees)
{
  // 0. Standardize names to 0..n-1 for this subproblem
  const auto& taxonomy = trees.back();
  auto all_leaves = taxonomy->getRoot()->getData().desIds;

  // index -> id
  vector<long> ids;
  // id -> index
  map<long,int> id_map;
  for(long id: all_leaves)
  {
    int i = ids.size();
    id_map[id] = i;
    ids.push_back(id);

    assert(id_map[ids[i]] == i);
    assert(ids[id_map[id]] == id);
  }
  auto remap = [&id_map](const set<long>& ids) {return remap_ids(ids,id_map);};
  vector<int> all_leaves_indices;
  for(int i=0;i<all_leaves.size();i++)
    all_leaves_indices.push_back(i);

  indices.resize(all_leaves.size());
  for(auto& i: indices)
    i=-1;
  
  // 1. Find splits in order of input trees
  vector<RSplit> splits;
  for(const auto& tree: trees)
  {
    auto root = tree->getRoot();
    const auto leafTaxa = root->getData().desIds;

    for(const auto& leaf: set_difference_as_set(leafTaxa, all_leaves))
      throw OTCError()<<"OTT Id "<<leaf<<" not in taxonomy!";
      
    const auto leafTaxaIndices = remap(leafTaxa);
    for(auto nd: iter_post_const(*tree))
    {
      const auto& descendants = remap(nd->getData().desIds);
      RSplit split{descendants, leafTaxaIndices};
      if (split.in.size()>1 and split.out.size())
      	splits.push_back(split);
    }
  }
  vector<const RSplit*> split_ptrs;
  for(const auto& split: splits)
    split_ptrs.push_back(&split);

  // 2. Add splits sequentially if they are consistent with previous splits.
  vector<const RSplit*> consistent;
  for(const auto& split: split_ptrs)
  {
    consistent.push_back(split);
    auto result = BUILD(all_leaves_indices, consistent);
    if (not result)
    {
      consistent.pop_back();
      LOG(DEBUG)<<"Reject: "<<*split<<"\n";
    }
    else
      LOG(DEBUG)<<"Keep:   "<<*split<<"\n";
  }

  // 3. Construct final tree and add names
  auto tree = BUILD(all_leaves_indices, consistent);
  for(auto nd: iter_pre(*tree))
    if (nd->isTip())
    {
      int index = nd->getOttId();
      nd->setOttId(ids[index]);
    }

  add_names(tree, taxonomy);
  return tree;
}

bool handleRequireOttIds(OTCLI & otCLI, const std::string & arg)
{
  if (arg == "true" or arg == "yes" or arg == "True" or arg == "Yes")
    otCLI.getParsingRules().requireOttIds = true;
  else if (arg == "false" or arg == "no" or arg == "False" or arg == "no")
    otCLI.getParsingRules().requireOttIds = false;
  else
    return false;
    
  return true;
}

bool handlePruneUnrecognizedTips(OTCLI & otCLI, const std::string & arg)
{
  if (arg == "true" or arg == "yes" or arg == "True" or arg == "Yes")
    otCLI.getParsingRules().pruneUnrecognizedInputTips = true;
  else if (arg == "false" or arg == "no" or arg == "False" or arg == "no")
    otCLI.getParsingRules().pruneUnrecognizedInputTips = false;
  else
    return false;
    
  return true;
}

bool synthesize_taxonomy = false;

bool handleSynthesizeTaxonomy(OTCLI &, const std::string & arg)
{
  if (arg == "true" or "yes" or "True" or "Yes")
    synthesize_taxonomy = true;
  else if (arg == "false" or "no" or "False" or "no")
    synthesize_taxonomy = false;
  else
    return false;
    
  return true;
}

unique_ptr<Tree_t> make_unresolved_tree(const vector<unique_ptr<Tree_t>>& trees, bool use_ids)
{
  std::unique_ptr<Tree_t> tree(new Tree_t());
  tree->createRoot();

  if (use_ids)
  {
    map<long,string> names;
    for(const auto& tree: trees)
      for(auto nd: iter_pre_const(*tree))
	if (nd->isTip())
	{
	  long id = nd->getOttId();
	  auto it = names.find(id);
	  if (it == names.end())
	    names[id] = nd->getName();
	}
  
    for(const auto& n: names)
    {
      auto node = tree->createChild(tree->getRoot());
      node->setOttId(n.first);
      node->setName(n.second);
    }
    clearAndfillDesIdSets(*tree);
  }
  else
  {
    set<string> names;
    for(const auto& tree: trees)
      for(auto nd: iter_pre_const(*tree))
	if (nd->isTip())
	  names.insert(nd->getName());

    for(const auto& n: names)
    {
      auto node = tree->createChild(tree->getRoot());
      node->setName(n);
    }
  }

  return tree;
}

string newick(const Tree_t &t)
{
  std::ostringstream s;
  writeTreeAsNewick(s, t);
  return s.str();
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-solve-subproblem",
                "Takes at series of tree files, with possibly mulitple trees per file.\n"
                "Files are concatenated and the combined list treated as a single subproblem.\n"
                "Trees should occur in order of priority, with the taxonomy last.\n",
		"subproblem.tre");

    otCLI.addFlag('o',
		  "Require OTT ids.  Defaults to true",
		  handleRequireOttIds,
		  true);

    otCLI.addFlag('p',
		  "Prune unrecognized tips.  Defaults to false",
		  handlePruneUnrecognizedTips,
		  true);

    otCLI.addFlag('T',
		  "Synthesize an unresolved taxonomy from all mentioned tips.  Defaults to false",
		  handleSynthesizeTaxonomy,
		  false);

    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};

    if (argc < 2)
	throw OTCError("No subproblem provided!");

    // I think multiple subproblem files are essentially concatenated.
    // Is it possible to read a single subproblem from cin?
    if (treeProcessingMain<Tree_t>(otCLI, argc, argv, get, nullptr, 1))
      std::exit(1);

    bool requireOttIds = otCLI.getParsingRules().requireOttIds;
    if (synthesize_taxonomy)
    {
      trees.push_back(make_unresolved_tree(trees,requireOttIds));
      LOG(DEBUG)<<"taxonomy = "<<newick(*trees.back())<<"\n";
    }
 
    // Add fake Ott Ids to tips and compute desIds
    if (not requireOttIds)
    {
      auto& taxonomy = trees.back();

      // 1. Compute mapping from name -> id
      long id = 1;
      map<string,long> name_to_id;
      for(auto nd: iter_pre(*taxonomy))
	if (nd->isTip())
	{
	  string name = nd->getName();
	  if (not name.size())
	    throw OTCError()<<"Taxonomy tip has no label or OTT id!";
	  auto it = name_to_id.find(name);
	  if (it != name_to_id.end())
	    throw OTCError()<<"Tip label '"<<name<<"' occurs twice in taxonomy!";
	  name_to_id[name] = id++;
	}

      // 2. Set ids
      for(auto& tree: trees)
	for(auto nd: iter_post(*tree))
	  if (nd->isTip())
	  {
	    string name = nd->getName();
	    if (not name.size())
	      throw OTCError()<<"Tip has no label or OTT id!";
	    auto it = name_to_id.find(name);
	    if (it == name_to_id.end())
	      throw OTCError()<<"Can't find label '"<<name<<"' in taxonomy!";
	    auto id = it->second;
	    nd->setOttId(id);
	    tree->getData().ottIdToNode[id] = nd;
	  }

      // 3. Compute DesIds.
      for(auto& tree: trees)
	clearAndfillDesIdSets(*tree);
    }

    auto tree = combine(trees);
    
    writeTreeAsNewick(std::cout, *tree);
    std::cout<<"\n";
}
