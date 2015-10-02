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
  set<long> in;
  set<long> out;
  set<long> all;
  RSplit() = default;
  RSplit(const set<long>& i, const set<long>& a)
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
long merge_components(long c1, long c2, map<long,long>& component, map<long,list<long>>& elements)
{
  if (elements[c2].size() > elements[c1].size())
    std::swap(c1,c2);

  for(long i:elements[c2])
    component[i] = c1;

  elements[c1].splice(elements[c1].begin(), elements[c2]);
  
  return c1;
}

bool empty_intersection(const set<long>& x, const set<long>& y)
{
  std::set<long>::const_iterator i = x.begin();
  std::set<long>::const_iterator j = y.begin();
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

/// Construct a tree with all the splits mentioned, and return a null pointer if this is not possible
unique_ptr<Tree_t> BUILD(const std::set<long>& tips, const vector<RSplit>& splits)
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
  map<long,long> component;      // element   -> component
  map<long,list<long>> elements; // component -> elements
  for(long i: tips)
  {
    component[i] = i;
    elements[i].push_back(i);
  }

  // 3. For each split, all the leaves in the include group must be in the same component
  for(const auto& split: splits)
  {
    long c1 = -1;
    for(long i: split.in)
    {
      long c2 = component[i];
      if (c1 != -1 and c1 != c2)
	merge_components(c1,c2,component,elements);
      c1 = component[i];
    }
  }

  // 4. If we can't subdivide the leaves in any way, then the splits are not consistent, so return failure
  long first = *tips.begin();
  if (elements[component[first]].size() == tips.size())
    return {};

  // 5. Create the set of tips in each connected component 
  map<long,set<long>> subtips;
  for(long c: tips)
  {
    if (c != component[c]) continue;
    
    set<long>& s = subtips[c];
    for(long l: elements[c])
      s.insert(l);
  }

  // 6. Determine the splits that are not satisfied yet and go into each component
  map<long,vector<RSplit>> subsplits;
  for(const auto& split: splits)
  {
    long first = *split.in.begin();
    long c = component[first];

    // if none of the exclude group are in the component, then the split is satisfied by the top-level partition.
    if (empty_intersection(split.out, subtips[c])) continue;

    subsplits[c].push_back(split);
  }
  
  // 7. Recursively solve the sub-problems of the partition components
  for(const auto& x: subtips)
  {
    auto subtree = BUILD(x.second, subsplits[x.first]);
    if (not subtree) return {};

    tree->addSubtree(tree->getRoot(), *subtree);
  }

  return tree;
}


/// Run the BUILD algorithm on the first n splits
unique_ptr<Tree_t> BUILD(const std::set<long>& tips, const vector<RSplit>& splits, int n)
{
  vector<RSplit> splits2;
  for(int i=0;i<n;i++)
    splits2.push_back(splits[i]);
  return BUILD(tips,splits2);
}

/// Copy node names from taxonomy to tree based on ott ids, and copy the root name also
void add_names(unique_ptr<Tree_t>& tree, const unique_ptr<Tree_t>& taxonomy)
{
  map<long,string> names;
  
  // 1. Determine the names for each ottid
  for(auto nd: iter_post(*taxonomy))
    if (nd->isTip())
      names[nd->getOttId()] = nd->getName();
  
  // 2. Set the names of nodes based on their ottid
  for(auto nd: iter_post(*tree))
    if (nd->hasOttId())
      nd->setName( names[nd->getOttId()] );
  
  // 3. Name the root node too
  string rootname = taxonomy->getRoot()->getName();
  tree->getRoot()->setName(rootname);
}

/// Get the list of splits, and add them one at a time if they are consistent with previous splits
unique_ptr<Tree_t> combine(const vector<unique_ptr<Tree_t>>& trees)
{
  // Standardize names to 0..n-1 for this subproblem
  const auto& taxonomy = trees.back();
  auto all_leaves = taxonomy->getRoot()->getData().desIds;
  
  // 1. Find splits in order of input trees
  vector<RSplit> splits;
  for(const auto& tree: trees)
  {
    auto root = tree->getRoot();
    const auto leafTaxa = root->getData().desIds;

    for(const auto& leaf: set_difference_as_set(leafTaxa, all_leaves))
      throw OTCError()<<"OTT Id "<<leaf<<" not in taxonomy!";
      
    for(auto nd: iter_post_const(*tree))
    {
      const auto& descendants = nd->getData().desIds;
      RSplit split{descendants, leafTaxa};
      if (split.in.size()>1 and split.out.size())
      	splits.push_back(split);
    }
  }

  // 2. Add splits sequentially if they are consistent with previous splits.
  vector<RSplit> consistent;
  for(const auto& split: splits)
  {
    consistent.push_back(split);
    auto result = BUILD(all_leaves, consistent);
    if (not result)
    {
      consistent.pop_back();
      LOG(DEBUG)<<"Reject: "<<split<<"\n";
    }
    else
      LOG(DEBUG)<<"Keep:   "<<split<<"\n";
  }

  // 3. Construct final tree and add names
  auto tree = BUILD(all_leaves, consistent);
  add_names(tree, taxonomy);
  return tree;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-solve-subproblem",
                "takes at series of tree files. Each is treated as a subproblem.\n",
		"subproblem.tre");

    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};

    // I don't think multiple subproblem files are essentially concatenated.
    // Is it possible to read a single subproblem from cin?
    treeProcessingMain<Tree_t>(otCLI, argc, argv, get, nullptr, 1);

    auto tree = combine(trees);
    
    writeTreeAsNewick(std::cout, *tree);
    std::cout<<"\n";
}
