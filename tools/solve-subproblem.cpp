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
  set<int> in;
  set<int> out;
  set<int> all;
  RSplit() = default;
  RSplit(const set<int>& i, const set<int>& a)
    :in(i),all(a)
  {
    std::set_difference(all.begin(),all.end(),in.begin(),in.end(),std::inserter(out,out.end()));
    assert(in.size() + out.size() == all.size());
  }
};

void merge(int c1, int c2, map<int,int>& component, map<int,list<int>>& elements)
{
  if (elements[c2].size() > elements[c1].size())
    std::swap(c1,c2);

  for(int i:elements[c2])
    component[i] = c1;

  elements[c1].splice(elements[c1].begin(), elements[c2]);
}

unique_ptr<Tree_t> BUILD(const std::set<int>& tips, const vector<RSplit>& splits)
{
  map<int,int> component;    // tip -> component
  map<int,list<int>> elements; // component -> elements
  for(int i: tips)
  {
    component[i] = i;
    elements[i].push_back(i);
  }

  for(const auto& split: splits)
  {
    int c1 = -1;
    for(int i: split.in)
    {
      int c2 = component[i];
      if (c1 != -1 and c1 != c2)
	merge(c1,c2,component,elements);
      c1 = component[i];
    }
  }

  int first = *tips.begin();
  if (elements[component[first]].size() == tips.size())
    std::cout<<"Failure: 1 component!\n";
  std::cout<<"Components:\n";
  for(const auto& l:elements)
  {
    if (l.second.empty()) continue;
    
    std::cout<<"   "<<l.second<<"\n";
  }
  return {};
}


unique_ptr<Tree_t> BUILD(const std::set<int>& tips, const vector<RSplit>& splits, int n)
{
  
  return {};
}

map<long,int> get_name_mapping(const vector<long>& ids)
{
  std::map<long,int> names;
  for(int i=0;i<ids.size();i++)
    names[ids[i]] = i;
  return names;
}

set<int> remap_names(const set<long>& ids, const map<long,int>& name)
{
  set<int> indices;
  for(auto id: ids)
  {
    auto i = name.find(id);
    if (i == name.end())
      throw OTCError()<<"Can't find ottid "<<id<<" in taxonomy";
    indices.insert(i->second);
  }
  return indices;
}

vector<long> get_vector(const set<long>& ids)
{
  vector<long> v;
  for(long id: ids)
    v.push_back(id);
  return v;
}

unique_ptr<Tree_t> merge(const vector<unique_ptr<Tree_t>>& trees)
{
  // Standardize names to 0..n-1 for this subproblem
  const auto& taxonomy = trees.back();
  auto leafset = taxonomy->getRoot()->getData().desIds;
  std::cout<<leafset<<std::endl;
  // i -> ids[i]
  vector<long> leaves = get_vector(leafset);
  // ids[i] -> i
  auto names = get_name_mapping(leaves);
  auto rename = [&names](const set<long>& ids) {return remap_names(ids,names);};
  auto all_leaves = rename(leafset);
  
  vector<RSplit> splits;
  for(const auto& tree: trees)
  {
    std::cout<<"Tree!\n";
    auto root = tree->getRoot();
    const auto leafTaxa = rename(root->getData().desIds);
    for(auto nd: iter_post_const(*tree))
    {
      const auto& descendants = rename(nd->getData().desIds);
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
