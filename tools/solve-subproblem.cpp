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

/// Merge components c1 and c2
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

  // 4. If we can't subdivide the leaves in any way, then the splits are not consistent.
  long first = *tips.begin();
  if (elements[component[first]].size() == tips.size())
  {
    std::cerr<<"Failure: 1 component!\n";
    return {};
  }

  // 5. Create the sets describing each component
  std::cerr<<"Components:\n";
  map<long,set<long>> subtips;
  for(long c: tips)
  {
    if (c != component[c]) continue;
    
    set<long>& s = subtips[c];
    for(long l: elements[c])
      s.insert(l);
    std::cerr<<"   "<<s<<"\n";
  }

  // 6. Determine the splits that are not satisfied yet and go into each component
  map<long,vector<RSplit>> subsplits;
  for(const auto& split: splits)
  {
    long first = *split.in.begin();
    long c = component[first];

    // if none of the exclude group are in the component, the the split is satisfied by the top-level partition.
    if (empty_intersection(split.out, subtips[c])) continue;

    subsplits[c].push_back(split);
  }
  
  // 7. Recursively solve the sub-problems of the partition components
  for(const auto& x: subtips)
  {
    int c = x.first;
    auto subtree = BUILD(x.second, subsplits[x.first]);
    if (not subtree) return {};

    tree->addSubtree(tree->getRoot(), *subtree);
  }

  return tree;
}


unique_ptr<Tree_t> BUILD(const std::set<long>& tips, const vector<RSplit>& splits, int n)
{
  vector<RSplit> splits2;
  for(int i=0;i<n;i++)
    splits2.push_back(splits[i]);
  return BUILD(tips,splits2);
}

string fixname(string name)
{
  for(char& c: name)
    if (c == ' ')
      c = '_';
  return name;
}

void add_names(unique_ptr<Tree_t>& tree, const unique_ptr<Tree_t>& taxonomy)
{
  map<long,string> names;
  
  // 8. Add writable names to the tree
  for(auto nd: iter_post(*taxonomy))
    if (nd->isTip())
      names[nd->getOttId()] = nd->getName();
  
  // 8. Add writable names to the tree
  for(auto nd: iter_post(*tree))
    if (nd->hasOttId())
    {
      string name = names[nd->getOttId()];
      nd->setName(fixname(name));
    }
  
  tree->getRoot()->setName(fixname(taxonomy->getRoot()->getName()));
}

unique_ptr<Tree_t> combine(const vector<unique_ptr<Tree_t>>& trees)
{
  // Standardize names to 0..n-1 for this subproblem
  const auto& taxonomy = trees.back();
  auto all_leaves = taxonomy->getRoot()->getData().desIds;
  std::cerr<<all_leaves<<std::endl;
  
  vector<RSplit> splits;
  for(const auto& tree: trees)
  {
    std::cerr<<"Tree!\n";
    auto root = tree->getRoot();
    const auto leafTaxa = root->getData().desIds;
    for(auto nd: iter_post_const(*tree))
    {
      const auto& descendants = nd->getData().desIds;
      RSplit split{descendants, leafTaxa};
      if (split.in.size()>1 and split.out.size())
      {
      	splits.push_back(split);
	std::cerr<<split.in<<" | "<<split.out<<"\n";
      }
    }
    BUILD(all_leaves, splits);
  }

  for(int n=2;n<=splits.size();)
  {
    unique_ptr<Tree_t> result = BUILD(all_leaves, splits, n);
    if (not result)
      splits.erase(splits.begin()+(n-1));
    else
      n++;
  }

  auto tree = BUILD(all_leaves, splits, splits.size());
  add_names(tree, taxonomy);
  return tree;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-subproblem-stats",
                "takes at series of tree files. Each is treated as a subproblem.\n" \
                "For each subproblem there is tab-separated report for:\n" \
                "    Subproblem name\n" \
                "    InSp = # of informative (nontrivial) splits\n" \
                "    LSS = size of the leaf label set\n" \
                "    ILSS = size of the set of labels included in at least one \"ingroup\"\n" \
                "    NT = The number of trees.\n" \
                "    TreeSummaryName = tree index or summary name\n" \
                "where the summary name can be Phylo-only or Total. \n" \
                "  \"Total\" summarizes info all trees in the file (including the taxonomy).\n" \
                "  \"Phylo-only\" former summarizes all of the phylogenetic inputs",
                "taxonomy.tre supertree/step_7_scratch/export-sub-temp/ott1000236.tre supertree/step_7_scratch/export-sub-temp/ott1000250.tre");

    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};

    // Want splits for each tip.
    // Not sure these trees are aligned to the taxonomy
    // They should have only terminals in desIds.
    treeProcessingMain<Tree_t>(otCLI, argc, argv, get, nullptr, 1);
    
    auto tree = combine(trees);
    
    tree->writeAsNewick(std::cout, true);
    std::cout<<";\n";
}
