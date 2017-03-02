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
void combine2(vector<unique_ptr<Tree_t>>& trees, bool verbose);
bool handleRequireOttIds(OTCLI & otCLI, const std::string & arg);
bool handlePruneUnrecognizedTips(OTCLI & otCLI, const std::string & arg);
bool handleRegraft(OTCLI&, const std::string & arg);
bool handleRootName(OTCLI&, const std::string & arg);


void combine2(vector<unique_ptr<Tree_t>>& trees, bool verbose) {
    assert(trees.size() == 2);
    auto& solution = *trees[0];
    auto& taxonomy = *trees[1];
    std::cerr << "Leaves:           solution = " << countLeaves(solution) << "   taxonomy = " << countLeaves(taxonomy) << std::endl;
    std::cerr << "Internal:         solution = " << n_internal(solution) << "   taxonomy = " << n_internal(taxonomy) << std::endl;
    std::cerr << "Internal splits:  solution = " << n_internal_out_degree_many(solution) << "   taxonomy = " << n_internal_out_degree_many(taxonomy) << std::endl;
    const auto out_degree_many1 = n_internal_out_degree_many(taxonomy);
    const auto n_internal_confirmed = n_internal_with_ott_id(solution);
    const auto n_internal_new = n_internal(solution) - n_internal_confirmed;
    // 1. First, remove nodes from the taxonomy that do not occur in the solution
    // 1a. Index solution nodes by OttId.
    map<long, Tree_t::node_type*> ott_to_sol;
    for (auto nd: iter_post(solution)){
        if (nd->hasOttId()){
            ott_to_sol[nd->getOttId()] = nd;
        }
    }
    map<long, Tree_t::node_type*> ott_to_tax;
    for (auto nd: iter_post(taxonomy)){
        if (nd->hasOttId()) {
            ott_to_tax[nd->getOttId()] = nd;
        } else {
            LOG(WARNING) << "  warning: node in taxonomy without an OTT ID.\n";
        }
    }
    for (auto nd: iter_post(solution)) {
        if (nd->isTip()) {
            if (not ott_to_tax.count(nd->getOttId())) {
                throw OTCError()<<"OttId "<<nd->getOttId()<<" not in taxonomy!";
            }
            //auto nd2 = ott_to_tax.at(nd->getOttId());
        }
    }
    // 1b. Find the subtree ancestral to the solution OttIds, and mark nodes monotypic in this subtree
    std::set<Tree_t::node_type*> ancestral;
    for (auto nd: iter_post(taxonomy)) {
        if (ott_to_sol.count(nd->getOttId())) {
            ancestral.insert(nd);
            auto a = nd->getParent();
            while (a) {
                ancestral.insert(a);
                a = a->getParent();
            }
        }
    }
    // 1c. Look at all ancestral nodes that are NOT monotypic
    //     Keep them if they OR one of their monotypic ancestors survives
    for (auto nd: all_nodes(taxonomy)) {
        if (not ancestral.count(nd)) {
            continue;
        }
        if (is_monotypic_in_set(nd, ancestral)) {
            continue;
        }
        Tree_t::node_type* nd1 = nullptr;
        if (ott_to_sol.count(nd->getOttId()) > 0) {
            nd1 = ott_to_sol.at(nd->getOttId());
        }
        vector<Tree_t::node_type*> nodes = {nd};
        auto anc = nd->getParent();
        while (not nd1 and anc and is_monotypic_in_set(anc,ancestral)) {
            nodes.push_back(anc);
            if (ott_to_sol.count(anc->getOttId()) > 0) {
                if (verbose){
                    LOG(INFO) << "Monotypic ancestor '" << anc->get_name() << "' in solution tree!";
                }
                nd1 = ott_to_sol.at(anc->getOttId());
            }
            anc = anc->getParent();
        }
        if (nd1) {
            if (not ott_to_sol.count(nd->getOttId())) {
                nd1 = bisect_branch_with_new_child(nd1);
                nd1->setOttId(nd->getOttId());
                nd1->setName(nd->get_name());
                ott_to_sol[nd1->getOttId()] = nd1;
            }
            assert(ott_to_sol.count(nd->getOttId()));
        } else {
            while(nodes.size()) {
                if (verbose) {
                    LOG(INFO)<<"Removing Id = '"<<nodes.back()->get_name()<<"' ("<<nodes.back()->getOttId()<<")"
                             <<"  children = "<<nodes.back()->getOutDegree()
                             <<"  ancestral children = "<<count_children_in_set(nodes.back(),ancestral);
                }
                // MTH this is where we should make note of which higher taxa do not make it into the solution.
                collapse_split_and_del_node(nodes.back());
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
    //map<const Tree_t::node_type*,Tree_t::node_type*> sol_to_tax;
    for (auto nd2: iter_post(taxonomy)) {
        if (ancestral.count(nd2) and not is_monotypic_in_set(nd2,ancestral)){
            assert(ott_to_sol.count(nd2->getOttId()));
        }
    }
    for (auto nd2: all_nodes(taxonomy)){
        if (ancestral.count(nd2) and ott_to_sol.count(nd2->getOttId())) {
            auto nd1 = ott_to_sol.at(nd2->getOttId());
            assert(nd1->getOttId() == nd2->getOttId());
            // Add the immediate ancestral nodes of nd2 to the solution tree, if they are monotypic
            while (nd2->getParent() and is_monotypic_in_set(nd2->getParent(),ancestral)
                   and not ott_to_sol.count(nd2->getParent()->getOttId())) {
                nd2 = nd2->getParent();
                assert(nd1->getParent());
                auto x = solution.createChild(nd1->getParent());
                nd1->detachThisNode();
                x->addChild(nd1);
                nd1 = x;
                nd1->setOttId(nd2->getOttId());
                nd1->setName(nd2->get_name());
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
            nd1->setName(nd2->get_name());
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
        throw OTCError() << "Supplied " << trees.size() << " trees for regrafting, should be 2 trees!";
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
