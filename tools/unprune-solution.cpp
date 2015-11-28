#include <algorithm>
#include <set>
#include <list>
#include <iterator>

#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
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
static bool regrafting = false;
static string rootName = "";


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

void remove_split(Tree_t::node_type* nd);
void show_children(Tree_t::node_type* nd, const set<Tree_t::node_type*>& ancestral);
int count_children(Tree_t::node_type* nd, const set<Tree_t::node_type*>& ancestral);
bool is_monotypic(Tree_t::node_type* nd, const set<Tree_t::node_type*>& ancestral);
Tree_t::node_type* insert_child(Tree_t::node_type* x);
long n_internal_with_ott_id(const Tree_t& T);
long n_internal(const Tree_t& T);
long n_internal_out_degree_1(const Tree_t& T);
long n_internal_out_degree_many(const Tree_t& T);
long n_nodes(const Tree_t& T);
std::vector<Tree_t::node_type*> all_nodes(Tree_t& T);
void combine2(vector<unique_ptr<Tree_t>>& trees, bool verbose);
bool handleRequireOttIds(OTCLI & otCLI, const std::string & arg);
bool handlePruneUnrecognizedTips(OTCLI & otCLI, const std::string & arg);
bool handleRegraft(OTCLI&, const std::string & arg);
bool handleRootName(OTCLI&, const std::string & arg);
unique_ptr<Tree_t> make_unresolved_tree(const vector<unique_ptr<Tree_t>>& trees, bool use_ids);

void remove_split(Tree_t::node_type* nd) {
    while (nd->getFirstChild()) {
        auto nd2 = nd->getFirstChild();
        nd2->detachThisNode();
        nd->addSibOnLeft(nd2);
    }
    nd->detachThisNode();
    delete nd;
}

void show_children(Tree_t::node_type* nd, const set<Tree_t::node_type*>& ancestral) {
    std::cerr<<nd->getName()<<" children: ";
    for (auto c = nd->getFirstChild(); c; c = c->getNextSib()) {
        if (ancestral.count(c)) {
            std::cerr << "'" << c->getName() << "' ";
        }
    }
    std::cerr << std::endl;
}

int count_children(Tree_t::node_type* nd, const set<Tree_t::node_type*>& ancestral) {
    int count = 0;
    for(auto c = nd->getFirstChild(); c; c = c->getNextSib()) {
        if (ancestral.count(c)) {
            count++;
        }
    }
    return count;
}

bool is_monotypic(Tree_t::node_type* nd, const set<Tree_t::node_type*>& ancestral) {
    return (count_children(nd,ancestral) == 1);
}

Tree_t::node_type* insert_child(Tree_t::node_type* x) {
    auto xc = new Tree_t::node_type(x);
    while (x->getFirstChild()) {
        auto nd = x->getFirstChild();
        nd->detachThisNode();
        xc->addChild(nd);
    }
    assert(not x->hasChildren());
    x->addChild(xc);
    return xc;
}

long n_internal_with_ott_id(const Tree_t& T) {
    long count = 0;
    for (auto nd: iter_post_const(T)) {
        if (not nd->isTip() and nd->hasOttId()) {
            count++;
        }
    }
    return count;
}

long n_internal(const Tree_t& T) {
    long count = 0;
    for (auto nd: iter_post_const(T)) {
        if (not nd->isTip()) {
            count++;
        }
    }
    return count;
}

long n_internal_out_degree_1(const Tree_t& T) {
    long count = 0;
    for (auto nd: iter_post_const(T)) {
        if (not nd->isTip() and nd->isOutDegreeOneNode()) {
            count++;
        }
    }
    return count;
}

long n_internal_out_degree_many(const Tree_t& T) {
    long count = 0;
    for(auto nd: iter_post_const(T)) {
        if (not nd->isTip() and not nd->isOutDegreeOneNode()) {
            count++;
        }
    }
    return count;
}

long n_nodes(const Tree_t& T) {
#pragma clang diagnostic ignored  "-Wunused-variable"
    long count = 0;
    for(auto nd: iter_post_const(T)){
        count++;
    }
    return count;
}

std::vector<Tree_t::node_type*> all_nodes(Tree_t& T) {
    std::vector<Tree_t::node_type*> nodes;
    for(auto nd: iter_post(T)) {
        nodes.push_back(nd);
    }
    return nodes;
}

void combine2(vector<unique_ptr<Tree_t>>& trees, bool verbose) {
    assert(trees.size() == 2);
    auto& solution = *trees[0];
    auto& taxonomy = *trees[1];
    std::cerr<<"Leaves:           solution = " << countLeaves(solution) <<"   taxonomy = " << countLeaves(taxonomy)<<std::endl;
    std::cerr<<"Internal:         solution = "<<n_internal(solution)<<"   taxonomy = "<<n_internal(taxonomy)<<std::endl;
    std::cerr<<"Internal splits:  solution = "<<n_internal_out_degree_many(solution)<<"   taxonomy = "<<n_internal_out_degree_many(taxonomy)<<std::endl;
    const auto out_degree_many1 = n_internal_out_degree_many(taxonomy);
    const auto n_internal_confirmed = n_internal_with_ott_id(solution);
    const auto n_internal_new = n_internal(solution) - n_internal_confirmed;
    // 1. First, remove nodes from the taxonomy that do not occur in the solution
    // 1a. Index solution nodes by OttId.
    map<long, Tree_t::node_type*> ott_to_sol;
    for(auto nd: iter_post(solution)){
        if (nd->hasOttId()){
            ott_to_sol[nd->getOttId()] = nd;
        }
    }
    map<long, Tree_t::node_type*> ott_to_tax;
    for(auto nd: iter_post(taxonomy)){
            if (nd->hasOttId()) {
                ott_to_tax[nd->getOttId()] = nd;
        }
    }
    for(auto nd: iter_post(solution)) {
        if (nd->isTip()) {
            if (not ott_to_tax.count(nd->getOttId())) {
                throw OTCError()<<"OttId "<<nd->getOttId()<<" not in taxonomy!";
            }
            auto nd2 = ott_to_tax.at(nd->getOttId());
        }
    }
    // 1b. Find the subtree ancestral to the solution OttIds, and mark nodes monotypic in this subtree
    std::set<Tree_t::node_type*> ancestral;
    for(auto nd: iter_post(taxonomy)) {
        if (ott_to_sol.count(nd->getOttId())) {
            ancestral.insert(nd);
            auto a = nd->getParent();
            while(a) {
                ancestral.insert(a);
                a = a->getParent();
            }
        }
    }
    // 1c. Look at all ancestral nodes that are NOT monotypic
    //     Keep them if they OR one of their monotypic ancestors survives
    for(auto nd: all_nodes(taxonomy)) {
        if (not ancestral.count(nd)) {
            continue;
        }
        if (is_monotypic(nd,ancestral)) {
            continue;
        }
        Tree_t::node_type* nd1 = nullptr;
        if (ott_to_sol.count(nd->getOttId()) > 0) {
            nd1 = ott_to_sol.at(nd->getOttId());
        }
        vector<Tree_t::node_type*> nodes = {nd};
        auto anc = nd->getParent();
        while(not nd1 and anc and is_monotypic(anc,ancestral)) {
            nodes.push_back(anc);
            if (ott_to_sol.count(anc->getOttId()) > 0) {
                if (verbose){
                    LOG(INFO) << "Monotypic ancestor '" << anc->getName() << "' in solution tree!";
                }
                nd1 = ott_to_sol.at(anc->getOttId());
            }
            anc = anc->getParent();
        }
        if (nd1) {
            if (not ott_to_sol.count(nd->getOttId())) {
                nd1 = insert_child(nd1);
                nd1->setOttId(nd->getOttId());
                nd1->setName(nd->getName());
                ott_to_sol[nd1->getOttId()] = nd1;
            }
            assert(ott_to_sol.count(nd->getOttId()));
        } else {
            while(nodes.size()) {
                if (verbose) {
                    LOG(INFO)<<"Removing Id = '"<<nodes.back()->getName()<<"' ("<<nodes.back()->getOttId()<<")"
                             <<"  children = "<<nodes.back()->getOutDegree()
                             <<"  ancestral children = "<<count_children(nodes.back(),ancestral);
                }
                remove_split(nodes.back());
                ancestral.erase(nodes.back());
                nodes.pop_back();
            }
        }
    }
    const auto out_degree_many2 = n_internal_out_degree_many(taxonomy);
    // CLAIM: Monotypic nodes can get removed from the tree, but monotypic nodes don't become polytypic,
    //        and polytypic nodes don't become monotypic.  Therefore we don't need to update the monotypic labels.
    
    // 2. Second, add nodes to the taxonomy from the solution
    // 2a. Map solution leaves to taxonomy leaves (walking up monotypic chimneys)
    map<const Tree_t::node_type*,Tree_t::node_type*> sol_to_tax;
    for(auto nd2: iter_post(taxonomy)) {
        if (ancestral.count(nd2) and not is_monotypic(nd2,ancestral)){
            assert(ott_to_sol.count(nd2->getOttId()));
        }
    }
    for(auto nd2: all_nodes(taxonomy)){
        if (ancestral.count(nd2) and ott_to_sol.count(nd2->getOttId())) {
            auto nd1 = ott_to_sol.at(nd2->getOttId());
            assert(nd1->getOttId() == nd2->getOttId());
            // Add nodes about nd2 to the solution tree, if they are monotypic
            while (nd2->getParent() and is_monotypic(nd2->getParent(),ancestral)
                   and not ott_to_sol.count(nd2->getParent()->getOttId())) {
                nd2 = nd2->getParent();
                assert(nd1->getParent());
                auto x = solution.createChild(nd1->getParent());
                nd1->detachThisNode();
                x->addChild(nd1);
                nd1 = x;
                nd1->setOttId(nd2->getOttId());
                nd1->setName(nd2->getName());
                ott_to_sol[nd1->getOttId()] = nd1;
            }
        }
    }
    for(auto nd2: iter_post(taxonomy)) {
        if (ancestral.count(nd2)) {
            assert(ott_to_sol.count(nd2->getOttId()));
        }
    }
    for(auto nd2: all_nodes(taxonomy)) {
        if (ancestral.count(nd2)) {
            auto id = nd2->getOttId();
            auto nd1 = ott_to_sol.at(id);
            nd1->setName(nd2->getName());
        } else {
            auto p2 = nd2->getParent();
            if (ancestral.count(p2)) {
                auto p1 = ott_to_sol.at(p2->getOttId());
                nd2->detachThisNode();
                p1->addChild(nd2);
            }
        }
    }
    // This is similar to, but different from, the number of non-monotypic nodes reject.
    // That is because the rejected nodes are marked as monotypic if they have no ANCESTRAL children.
    std::cerr<<"Taxonomy splits: #rejected  by phylo inputs = "<<out_degree_many1 - out_degree_many2<<std::endl;
    std::cerr<<"Solution splits: #in taxonomy            = "<<n_internal_confirmed<<std::endl;
    std::cerr<<"Solution splits: #from phylo inputs only = "<<n_internal_new<<std::endl;
    const auto out_degree_many3 = n_internal_out_degree_many(solution);
    std::cerr<<"Unpruned splits: #added by phylo inputs = "<<out_degree_many3 - out_degree_many2<<std::endl;
    std::cerr<<"Unpruned splits: total = "<<out_degree_many3<<std::endl;
}

bool handleRequireOttIds(OTCLI & otCLI, const std::string & arg) {
    otCLI.getParsingRules().setOttIds = get_bool(arg,"-o: ");
    return true;
}

bool handlePruneUnrecognizedTips(OTCLI & otCLI, const std::string & arg) {
    otCLI.getParsingRules().pruneUnrecognizedInputTips = get_bool(arg,"-p: ");
    return true;
}

bool handleRegraft(OTCLI&, const std::string & arg) {
    if (arg.size()) {
        throw OTCError()<<"-r does not take an argument.";
    }
    regrafting = true;
    return true;
}


bool handleRootName(OTCLI&, const std::string & arg) {
    rootName = arg;
    return true;
}


/// Create an unresolved taxonomy out of all the input trees.
unique_ptr<Tree_t> make_unresolved_tree(const vector<unique_ptr<Tree_t>>& trees, bool use_ids) {
    std::unique_ptr<Tree_t> retTree(new Tree_t());
    retTree->createRoot();
    if (use_ids) {
        map<long,string> names;
        for(const auto& tree: trees){
            for(auto nd: iter_pre_const(*tree)){
                if (nd->isTip()) {
                    long id = nd->getOttId();
                    auto it = names.find(id);
                    if (it == names.end()) {
                        names[id] = nd->getName();
                    }
                }
            }
        }
        for(const auto& n: names) {
            auto node = retTree->createChild(retTree->getRoot());
            node->setOttId(n.first);
            node->setName(n.second);
        }
    } else {
        set<string> names;
        for(const auto& tree: trees) {
            for(auto nd: iter_pre_const(*tree)){
                if (nd->isTip()) {
                    names.insert(nd->getName());
                }
            }
        }
        for(const auto& n: names) {
            auto node = retTree->createChild(retTree->getRoot());
            node->setName(n);
        }
    }
    return retTree;
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-solve-subproblem",
                "Takes a series of tree files.\n"
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
    otCLI.addFlag('n',
                  "Rename the root to this name",
                  handleRootName,
                  true);
    otCLI.addFlag('r',
                  "Regrafting pruned leaves: assume two trees, a solution tree and\n"
                  "a taxonomy.  Determine conflicting clades in the taxonomy based on\n"
                  "which OTT Ids occur in the solution tree.",
                  handleRegraft,
                  false);
    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};
    if (argc < 2){
        throw OTCError("No subproblem provided!");
    }
    // I think multiple subproblem files are essentially concatenated.
    // Is it possible to read a single subproblem from cin?
    if (treeProcessingMain<Tree_t>(otCLI, argc, argv, get, nullptr, 1)){
        std::exit(1);
    }

    if (trees.empty()) {
        throw OTCError("No trees loaded!");
    }
    bool setOttIds = otCLI.getParsingRules().setOttIds;
    // Add fake Ott Ids to tips
    if (not setOttIds) {
        auto name_to_id = createIdsFromNames(*trees.back());
        for(auto& tree: trees) {
            setIdsFromNames(*tree, name_to_id);
        }
    }
    if (trees.size() != 2) {
        throw OTCError()<<"Supplied "<<trees.size()<<" trees for regrafting, should be 2 trees!";
    }
    combine2(trees, otCLI.verbose);
    unique_ptr<Tree_t> tree = std::move(trees[0]);
    if (not rootName.empty()){
        tree->getRoot()->setName(rootName);
    }
    writeTreeAsNewick(std::cout, *tree);
    std::cout<<"\n";
    return 0;
}
