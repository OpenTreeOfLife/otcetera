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
    std::set_difference(all.begin(),all.end(),in.begin(),in.end(),std::inserter(out,out.end()));
    assert(in.size() + out.size() == all.size());
  }
};

void merge(long c1, long c2, map<long,long>& component, map<long,list<long>>& elements)
{
  if (elements[c2].size() > elements[c1].size())
    std::swap(c1,c2);

  for(long i:elements[c2])
    component[i] = c1;

  elements[c1].splice(elements[c1].begin(), elements[c2]);
}

unique_ptr<Tree_t> BUILD(const std::set<long>& tips, const vector<RSplit>& splits)
{
  std::unique_ptr<Tree_t> tree(new Tree_t());
  tree->createRoot();

  if (tips.size() == 1)
  {
    auto Node1 = tree->createChild(tree->getRoot());
    Node1->setOttId(*tips.begin());
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

  map<long,long> component;    // tip -> component
  map<long,list<long>> elements; // component -> elements
  for(long i: tips)
  {
    component[i] = i;
    elements[i].push_back(i);
  }

  for(const auto& split: splits)
  {
    long c1 = -1;
    for(long i: split.in)
    {
      long c2 = component[i];
      if (c1 != -1 and c1 != c2)
	merge(c1,c2,component,elements);
      c1 = component[i];
    }
  }

  long first = *tips.begin();
  if (elements[component[first]].size() == tips.size())
  {
    std::cout<<"Failure: 1 component!\n";
    return {};
  }

  std::cout<<"Components:\n";
  map<long,set<long>> subtips;
  for(long c: tips)
  {
    if (c != component[c]) continue;
    
    set<long>& s = subtips[c];
    for(long l: elements[c])
      s.insert(l);
    std::cout<<"   "<<s<<"\n";
  }

  map<long,vector<RSplit>> subsplits;
  for(const auto& split: splits)
  {
    long first = *split.in.begin();
    long c = component[first];
    
  }
  
  return tree;
}


unique_ptr<Tree_t> BUILD(const std::set<long>& tips, const vector<RSplit>& splits, int n)
{
  
  return {};
}

unique_ptr<Tree_t> merge(const vector<unique_ptr<Tree_t>>& trees)
{
  // Standardize names to 0..n-1 for this subproblem
  const auto& taxonomy = trees.back();
  auto all_leaves = taxonomy->getRoot()->getData().desIds;
  std::cout<<all_leaves<<std::endl;
  
  vector<RSplit> splits;
  for(const auto& tree: trees)
  {
    std::cout<<"Tree!\n";
    auto root = tree->getRoot();
    const auto leafTaxa = root->getData().desIds;
    for(auto nd: iter_post_const(*tree))
    {
      const auto& descendants = nd->getData().desIds;
      RSplit split{descendants, leafTaxa};
      if (split.in.size()>1 and split.out.size())
      {
      	splits.push_back(split);
	std::cout<<split.in<<" | "<<split.out<<"\n";
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

  return BUILD(all_leaves, splits, splits.size());
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
    
    merge(trees);
}
