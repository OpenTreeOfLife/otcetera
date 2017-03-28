#include <queue>
#include "otc/greedy_forest.h"
#include "otc/node_embedding.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/debug.h"
#include "otc/supertree_context.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/error.h"
namespace otc {

template<typename T, typename U>
bool NodeEmbedding<T, U>::debug_node_embeddings(const char * tag, 
                                             bool is_contested,
                                             const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const {
    for (const auto  & t2exit : edgeBelowEmbeddings) {
        const auto treeIndex = t2exit.first;
        const auto & exitPaths =  t2exit.second;
        for (const auto epref : exitPaths) {
            const auto & sNERef = sn2ne.at(epref->scaffold_des);
            assert(contains(sNERef.edgeBelowEmbeddings, treeIndex));
            assert(contains(sNERef.edgeBelowEmbeddings.at(treeIndex), epref));
            for (auto aNode : iter_anc(*(epref->scaffold_des))) {
                if (aNode == epref->scaffold_anc) {
                    break;
                }
                const auto & saNERef = sn2ne.at(epref->scaffold_des);
                assert(contains(saNERef.edgeBelowEmbeddings, treeIndex));
                assert(contains(saNERef.edgeBelowEmbeddings.at(treeIndex), epref));
            }
        }
    }
    for (const auto & t2l : loopEmbeddings) {
        auto treeIndex = t2l.first;
        const auto childExitForThisTree = get_all_child_exit_paths_for_tree(embeddedNode,
                                                                      treeIndex,
                                                                      sn2ne);
        auto lnd2par = get_looped_phylo_node_to_par(treeIndex);
        auto end2par = get_exit_phylo_node_to_par(treeIndex);
        // all of the loop phyloParents should be in a child of another
        //  loop or a child of an exit path.
        // If the node is uncontested, then 
        U * root = nullptr;
        std::set<U *> parentsOfExits;
        for (const auto & c2p : lnd2par) {
            auto p = c2p.second;
            if (p != root && !contains(lnd2par, p)) {
                if (contains(end2par, p)) {
                    parentsOfExits.insert(end2par.at(p));
                    assert(root == nullptr || !parentsOfExits.empty());
                    assert(is_contested || parentsOfExits.size() < 2);
                } else {
                    if (p->get_parent() != nullptr) {
                        auto gp = p->get_parent();
                        LOG(ERROR)  << " Failing on scaffold ott" << embeddedNode->get_ott_id() << " treeIndex " << treeIndex;
                        LOG(ERROR)  << " phylo parent node (@" << reinterpret_cast<long>(p) << ") = " << get_designator(*p) << ' ' << reinterpret_cast<long>(p) << " is not in loop or exit edges";
                        write_newick(std::cerr, p);
                        std::cerr << std::endl;
                        LOG(ERROR)  << " the phylo child was (@" << reinterpret_cast<long>(c2p.first) << ") = " << get_designator(*c2p.first);
                        write_newick(std::cerr, c2p.first);
                        LOG(ERROR)  << " phylo grandparent (@" << reinterpret_cast<long>(gp) << ") = " << get_designator(*gp) << ' ' << reinterpret_cast<long>(gp);
                        write_newick(std::cerr, gp);
                        std::cerr << std::endl;
                        for (auto anc : iter_anc(*embeddedNode)) {
                            auto pathPairings = sn2ne.at(anc).refers_to_node(treeIndex, gp);
                            if (!pathPairings.empty()) {
                                LOG(ERROR) << " ancestor ott" << anc->get_ott_id() << " refers to the parent";
                            }
                        }
                        debugPrint(*embeddedNode, treeIndex, sn2ne);
                        LOG(ERROR) << "Failure context is " << tag;
                        assert(false);
                    }
                    assert(parentsOfExits.empty());
                }
                root = p;
            }
        }
        if (parentsOfExits.empty()) {
            // no loops, the tree is just polytomy for this subproblem
            for (auto pathPtr : childExitForThisTree) {
                const auto & rids = pathPtr->get_ott_id_set();
                assert(rids.size() == 1);
                assert(*rids.begin() != LONG_MAX);
            }
        } else {
            for (auto pp : childExitForThisTree) {
                if (pp->phylo_parent != root) {
                    parentsOfExits.insert(pp->phylo_parent);
                    assert(is_contested || parentsOfExits.size() < 2);
                }
                assert(!contains(lnd2par, pp->phylo_child));
                if (contains(end2par, pp->phylo_child)) {
                    assert(!contains(lnd2par, pp->phylo_parent));
                } else {
                    assert(contains(lnd2par, pp->phylo_parent));
                }
                const auto & rids = pp->get_ott_id_set();
                assert(rids.size() == 1);
                assert(rids.size() == 1);
                assert(*rids.begin() != LONG_MAX);
            }
        }
    }
    return true;
}

template<typename T, typename U>
void NodeEmbedding<T, U>::set_ott_id_for_exit_embeddings(
                    T * newScaffDes,
                    OttId ottId,
                    std::map<const T *, NodeEmbedding<T, U> > & n2ne) {
    for (auto treeInd2eout : edgeBelowEmbeddings) {
        assert(treeInd2eout.second.size() < 2);
        for (auto eout : treeInd2eout.second) {
            LOG(DEBUG) << "for tree " << treeInd2eout.first << " set_ott_id(" << ottId<< ')';
            eout->scaffold_des = newScaffDes;
            eout->set_ott_id_set(ottId, n2ne);
        }
    }
}

template<typename T, typename U>
void NodeEmbedding<T, U>::merge_exit_embeddings_if_multiple() {
    std::map<std::size_t, std::set<PathPairPtr> > toCull;
    for (auto treeInd2eout : edgeBelowEmbeddings) {
        if (treeInd2eout.second.size() > 1) {
            OTC_UNREACHABLE; // now resolving earlier...
            // If an input tree has a polytomy with members of a taxon as well as its "outgroup" taxa,
            //  then the polytomy does not contest the monophyly.
            //  this function will be called after resolution of the polytomy. 
            //  TODO: We might want to add a node for the resolved node should be added to the source tree to reflect
            //      that its ambiguity has been resolved in one particular way.
            //  The crucial thing for the continuation of the export/supertree operation
            //      is the set of exit nodes from embeddedNode be reduced to 1 exit node
            U * resolvedPolyParent = nullptr;
            for (auto eout : treeInd2eout.second) {
                if (resolvedPolyParent == nullptr) {
                    resolvedPolyParent = eout->phylo_parent;
                } else {
                    assert(resolvedPolyParent == eout->phylo_parent);
                    toCull[treeInd2eout.first].insert(eout); //TMP just keeping the first node
                }
            }
        }
    }
    for (const auto & toCullEl : toCull) {
        for (auto doomedPath : toCullEl.second) {
            remove_ref_to_exit_path(toCullEl.first, doomedPath);
        }
    }
}

template<typename T, typename U>
void NodeEmbedding<T, U>::resolve_parent_in_favor_of_this_node(
                T & scaffold_node,
                std::size_t treeIndex,
                SupertreeContextWithSplits & sc) {
    auto exitForThisTreeIt = edgeBelowEmbeddings.find(treeIndex);
    if (exitForThisTreeIt == edgeBelowEmbeddings.end()) {
        return;
    }
    PathPairSet & exitSetForThisTree{exitForThisTreeIt->second};
    const auto origNumExits = exitSetForThisTree.size();
    if (origNumExits < 2) {
        return;
    }
    const TreeMappedWithSplits * cTreePtr = sc.trees_by_index.at(treeIndex);
    // this is the only way that we modify the embedded tree. So for now, we'll
    //      tolerate the const cast, which is in very poor form...
    // @TMP @TODO
    TreeMappedWithSplits * treePtr = const_cast<TreeMappedWithSplits *>(cTreePtr);
    LOG(INFO) << "Resolving treeIndex = " << treeIndex << " name = " << treePtr->get_name() << " for OTT " << scaffold_node.get_ott_id() << '\n';
    const std::map<U *, U *> phyloNode2PhyloPar = get_exit_phylo_node_to_par(treeIndex);
    // resolve the phylo tree
    U * phPar = nullptr;
    for (auto p2p : phyloNode2PhyloPar) {
        if (phPar == nullptr) {
            phPar = p2p.second;
        } else {
            assert(phPar == p2p.second); // there should only be one par.
        }
    }
    assert(phPar != nullptr);
    U * insertedNodePtr = treePtr->create_node(phPar);
    for (auto p2p : phyloNode2PhyloPar) {
        U * phChild = p2p.first;
        phChild->detach_this_node();
        insertedNodePtr->add_child(phChild);
        const OttIdSet & cd = phChild->get_data().des_ids;
        insertedNodePtr->get_data().des_ids.insert(cd.begin(), cd.end());
    }
    T * scaffoldAncestor = nullptr;
    for (auto epp : exitSetForThisTree) {
        if (scaffoldAncestor == nullptr) {
            scaffoldAncestor = epp->scaffold_anc;
        } else {
            assert(scaffoldAncestor == epp->scaffold_anc); // to be a resolution case all exits must have the same parent
        }
    }
    // create the new node pairing and path pairing objects
    sc.node_pairings_from_resolve.emplace_back(NodePairingWithSplits(&scaffold_node, insertedNodePtr));
    NodePairingWithSplits & newNodePairing{*sc.node_pairings_from_resolve.rbegin()};
    nodeEmbeddings[treeIndex].insert(&newNodePairing);
    sc.path_pairings_from_resolve.emplace_back(PathPairingWithSplits(scaffoldAncestor, phPar, newNodePairing));
    PathPairingWithSplits & newPathPairing{*sc.path_pairings_from_resolve.rbegin()};
    // fix every exit path to treat scaffold_node as the scaffold_anc node.
    // If the scaffold_des is the scaffold_node, then this exit is becoming a loop...
    PathPairSet toMoveToLoops;
    std::map<const T *, NodeEmbedding<T, U> > & sn2ne = sc.scaffold_to_node_embedding;
    for (auto epp : exitSetForThisTree) {
        assert(epp->phylo_parent == phPar);
        epp->phylo_parent = insertedNodePtr;
        assert(epp->scaffold_anc != &scaffold_node);
        for (auto n : iter_anc(scaffold_node)) {
            if (n == epp->scaffold_anc) {
                break;
            }
            sn2ne.at(n).remove_ref_to_exit_path(treeIndex, epp);
        }
        epp->scaffold_anc = &scaffold_node;
        if (epp->scaffold_des == &scaffold_node) {
            toMoveToLoops.insert(epp);
        }
    }
    // all of the former exit paths for this tree are invalid...
    exitSetForThisTree.clear();
    // add the new loops...
    if (!toMoveToLoops.empty()) {
        auto laIt = loopEmbeddings.find(treeIndex);
        if (laIt == loopEmbeddings.end()) {
            loopEmbeddings[treeIndex] = toMoveToLoops;
        } else {
            laIt->second.insert(begin(toMoveToLoops), end(toMoveToLoops));
        }
    }
    // insert the new (and only) exit path for this node and its ancestors...
    edgeBelowEmbeddings[treeIndex].insert(&newPathPairing);
    for (auto n : iter_anc(scaffold_node)) {
        if (n == scaffoldAncestor) {
            break;
        }
        sn2ne.at(n).edgeBelowEmbeddings[treeIndex].insert(&newPathPairing);
    }
    assert(edgeBelowEmbeddings[treeIndex].size() == 1);
}



// Returns all loop paths for nd and all edgeBelowEmbeddings of its children
template<typename T, typename U>
std::vector<const PathPairing<T, U> *>
NodeEmbedding<T, U>::get_all_incoming_path_pairs(const std::map<const T *, NodeEmbedding<T, U> > & eForNd,
                                                 std::size_t treeIndex) const {
    const T *nd = embeddedNode;
    std::vector<const PathPairingWithSplits *> r;
    const auto lait = loopEmbeddings.find(treeIndex);
    if (lait != loopEmbeddings.end()) {
        for (const auto pp : lait->second) {
            r.push_back(pp);
        }
    }
    for (auto c : iter_child_const(*nd)) {
        //LOG(DEBUG) << "    get_all_incoming_path_pairs c = " << get_designator(*c);
        const auto cembed = eForNd.find(c);
        if (cembed == eForNd.end()) {
            //LOG(DEBUG) << "     No embedding found";
            continue;
        }
        const auto & emb = cembed->second;
        const auto ceait = emb.edgeBelowEmbeddings.find(treeIndex);
        if (ceait != emb.edgeBelowEmbeddings.end()) {
            for (const auto & e : ceait->second) {
                r.push_back(e);
            }
        } else {
                //LOG(DEBUG) << "     No embedding found for this tree";
        }
    }
    return r;
}

template<typename T, typename U>
bool PathPairing<T, U>::update_ott_id_set_no_traversal(const OttIdSet & oldEls, const OttIdSet & newEls) {
    if (false && debugging_output_enabled) {
        LOG(DEBUG) << "  update_ott_id_set_no_traversal for " << reinterpret_cast<long>(this) << " in ";
        db_write_ott_id_set("curr_child_ott_id_set", curr_child_ott_id_set);
        db_write_ott_id_set("oldEls", oldEls);
        db_write_ott_id_set("newEls", newEls);
    }
    auto i = set_intersection_as_set(oldEls, curr_child_ott_id_set);
    if (i.size() < oldEls.size()) {
        return false;
    }
    if (!i.empty()) {
        for (auto o : i) {
            curr_child_ott_id_set.erase(o);
        }
    }
    curr_child_ott_id_set.insert(begin(newEls), end(newEls));
    if (false && debugging_output_enabled) {
        LOG(DEBUG) << "  update_ott_id_set_no_traversal for " << reinterpret_cast<long>(this);
        db_write_ott_id_set("update_ott_id_set_no_traversal exit ", curr_child_ott_id_set);
    }
    return true;
}

template<typename T, typename U>
inline const OttIdSet & NodeEmbedding<T, U>::get_relevant_des_ids_from_path(const PathPairing<T, U> & path) {
    return path.get_ott_id_set();
}
template<typename T, typename U>
OttIdSet NodeEmbedding<T, U>::get_relevant_des_ids_from_path_pair_set(const PathPairSet & pps) {
    OttIdSet relevantIds;
    for (auto path : pps) {
        const auto & cdi = get_relevant_des_ids_from_path(*path);
        relevantIds.insert(begin(cdi), end(cdi));
    }
    return relevantIds;
}

template<typename T, typename U>
void NodeEmbedding<T, U>::collapse_source_edge(const T * , //phylo_parent,
                                                 PathPairing<T, U> * ) { //path
    assert("not implemented"[0] == 'f');; // until we check for "high ranks preserve contested monophyly optimization"
}

template<typename T, typename U>
void NodeEmbedding<T, U>::collapse_source_edge_to_force_one_entry(T & ,
                                                             PathPairSet & pps,
                                                             std::size_t treeIndex,
                                                             SupertreeContextWithSplits & sc) {
    assert("not implemented"[0] == 'f');; // until we check for "high ranks preserve contested monophyly optimization"
    if (pps.size() < 2) {
        return;
    }
    auto relevantIds = get_relevant_des_ids(sc.scaffold_to_node_embedding, treeIndex);
    PathPairing<T, U> * firstPairing = *pps.begin();
    const T * onePhyloPar = firstPairing->phylo_parent;
    const T * phyloMrca = search_anc_for_mrca_of_des_ids(onePhyloPar, relevantIds);
    std::set<const T *> prevCollapsed; 
    prevCollapsed.insert(phyloMrca); // we don't actually collapse this edge, we just add it to the set so we don't collapse it below....
    for (auto path : pps) {
        const auto pp = path->phylo_parent;
        if (!contains(prevCollapsed, pp)) {
            collapse_source_edge(pp, path);
            prevCollapsed.insert(pp);
        }
    }
}
template<typename T, typename U>
void NodeEmbedding<T, U>::resolve_given_contested_monophyly(T & scaffold_node,
                                                         SupertreeContextWithSplits & sc) {
    for (std::size_t treeIndex = 0 ; treeIndex < sc.num_trees; ++treeIndex) {
        const auto ebaIt = edgeBelowEmbeddings.find(treeIndex);
        if (ebaIt == edgeBelowEmbeddings.end()) {
            continue;
        }
        PathPairSet & pps = ebaIt->second;
        collapse_source_edge_to_force_one_entry(scaffold_node, pps, treeIndex, sc);
    }
    resolve_given_uncontested_monophyly(scaffold_node, sc);
}

template<typename T, typename U>
std::set<PathPairing<T, U> *> NodeEmbedding<T, U>::get_all_child_exit_paths(
                const T & scaffold_node,
                const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const {
    std::set<PathPairing<T, U> *> r;
    for (auto c : iter_child_const(scaffold_node)) {
        const auto & thr = sn2ne.at(c);
        for (auto te : thr.edgeBelowEmbeddings) {
            r.insert(begin(te.second), end(te.second));
        }
    }
    return r;
}

enum DOTFileStep {
    INFORMATIVE_SPLIT,
    TRIVIAL_SPLIT,
    NONEMBEDDED_SPLIT
};
const std::string getForestDOTFilename(const std::string & prefix,
                                       const DOTFileStep step,
                                       std::size_t treeIndex,
                                       OttId groupIndex);
const std::string getForestDOTFilename(const std::string & prefix,
                                       const DOTFileStep step,
                                       std::size_t treeIndex,
                                       OttId groupIndex) {
    std::string base = prefix;
    if (step == INFORMATIVE_SPLIT) {
        base += "AfterInfTree";
    } else if (step == TRIVIAL_SPLIT) {
        base += "AfterTrivTree";
    } else if (step == NONEMBEDDED_SPLIT) {
        base += "AfterXChildSplit";
    }
    if (step == INFORMATIVE_SPLIT || step == TRIVIAL_SPLIT || step == NONEMBEDDED_SPLIT) {
        base += std::to_string(treeIndex);
        base += "Split";
        base += std::to_string(groupIndex);
    }
    base += ".dot";
    return base;
}

template<typename T, typename U>
std::set<PathPairing<T, U> *> NodeEmbedding<T, U>::get_all_child_exit_paths_for_tree(
                const T & scaffold_node,
                std::size_t treeIndex,
                const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const {
    std::set<PathPairing<T, U> *> r;
    for (auto c : iter_child_const(scaffold_node)) {
        const auto & thr = sn2ne.at(c);
        const auto & tebeIt = thr.edgeBelowEmbeddings.find(treeIndex);
        if (tebeIt != thr.edgeBelowEmbeddings.end()) {
            r.insert(begin(tebeIt->second), end(tebeIt->second));
        }
    }
    return r;
}

template<typename T>
void debugPrintNd2Par(const char * p, const T & m) {
    std::cerr << "debugPrintNd2Par " << p << ": ";
    for (auto mp : m) {
        std::cerr << "  " << reinterpret_cast<long>(mp.first) << " -> " << reinterpret_cast<long>(mp.second) << '\n';
    }
}
// Only resolves cases in which there is a deeper polytomy that is compatible with 
//  the taxon - resolved to make the taxon monophyletic...
// This allows export of deeper subproblems to be performed.
// returns the name (not the full path) of the tree file written or an empty string
template<typename T, typename U>
std::string NodeEmbedding<T, U>::export_subproblem_and_resolve(
                            T & scaffold_node,
                            const std::string & exportDir,
                            std::ostream * exportStream,
                            SupertreeContextWithSplits & sc) {
    const std::map<const T *, NodeEmbedding<T, U> > & sn2ne = sc.scaffold_to_node_embedding;
    //debug_node_embeddings("top of export", false, sn2ne);
    //debugPrint(scaffold_node, 215, sn2ne);
    const OttIdSet EMPTY_SET;
    const auto scaffOTTId = scaffold_node.get_ott_id();
    std::ofstream treeFileStream;
    std::ofstream provFileStream;
    std::ostream * treeExpStream = exportStream;
    std::ostream * provExpStream = exportStream;
    std::string retStr;
    if (exportStream == nullptr) {
        std::string outFilename = exportDir;
        outFilename.append("/ott");
        outFilename += std::to_string(scaffOTTId);
        outFilename += ".tre";
        retStr.append("ott");
        retStr += std::to_string(scaffOTTId);
        retStr += ".tre";
        std::string provOutFilename = exportDir;
        provOutFilename.append("/ott");
        provOutFilename += std::to_string(scaffOTTId);
        provOutFilename += "-tree-names.txt";
        treeFileStream.open(outFilename);
        if (!treeFileStream.good()) {
            std::string m = "Could not open \"";
            m += outFilename;
            m += '\"';
            throw OTCError(m);
        }
        treeExpStream = &treeFileStream;
        provFileStream.open(provOutFilename);
        if (!provFileStream.good()) {
            std::string m = "Could not open \"";
            m += provOutFilename;
            m += '\"';
            throw OTCError(m);
        }
        provExpStream = &provFileStream;
    } else {
        *exportStream << "ott" << scaffOTTId << '\n';
    }
    
    //TMP this could be done a lot more efficiently. Writing the trees through the GBF should
    // exercise some of the code for the scaffolded supertree operation. So this should help
    //  us find bugs in that code...
    
    OttIdSet totalLeafSet;
    for (std::size_t treeIndex = 0 ; treeIndex < sc.num_trees; ++treeIndex) {
        const auto * treePtr = sc.trees_by_index.at(treeIndex);
        assert(treePtr != nullptr);
        //auto prelnd2par = get_looped_phylo_node_to_par(treeIndex);
        //debugPrintNd2Par("preloops", prelnd2par);
        //auto preend2par = get_exit_phylo_node_to_par(treeIndex);
        //debugPrintNd2Par("preexits", preend2par);
        resolve_parent_in_favor_of_this_node(scaffold_node, treeIndex, sc);
        const auto childExitForThisTree = get_all_child_exit_paths_for_tree(scaffold_node,
                                                                      treeIndex,
                                                                      sn2ne);
        const auto tempDBIt = edgeBelowEmbeddings.find(treeIndex);
        assert(tempDBIt == edgeBelowEmbeddings.end() || tempDBIt->second.size() < 2);
        auto laIt = loopEmbeddings.find(treeIndex);
        if (laIt == loopEmbeddings.end() && childExitForThisTree.empty()) {
            continue;
        }
        const auto firstLaIt = laIt;
        LOG(INFO) << "     exporting ott" << scaffOTTId << "  treeIndex = " << treeIndex << " name = " << treePtr->get_name();
        const OttIdSet relevantIds = get_relevant_des_ids(sn2ne, treeIndex);
        totalLeafSet.insert(begin(relevantIds), end(relevantIds));
        auto lnd2par = get_looped_phylo_node_to_par(treeIndex);
        //debugPrintNd2Par("loops", lnd2par);
        const auto shouldHaveBeenPrunedNd2Par = get_un_embedded_phylo_node_to_par(treeIndex);
        std::map<NodeWithSplits *, OttId> nd2id;
        OttIdSet shouldHaveBeenPrunedIds;
        for (auto unembeddedPair : shouldHaveBeenPrunedNd2Par) {
            auto uDes = unembeddedPair.first;
            auto uPar = unembeddedPair.second;
            assert(!contains(lnd2par, uDes));
            lnd2par[uDes] = uPar;
            auto shbpId = uDes->get_ott_id();
            nd2id[uDes] = shbpId;
            shouldHaveBeenPrunedIds.insert(shbpId);
        }
        
        auto end2par = get_exit_phylo_node_to_par(treeIndex);
        //debugPrintNd2Par("exits", end2par);
        U * root = nullptr;
        std::set<U *> parentsOfExits;
        for (const auto & c2p : lnd2par) {
            auto p = c2p.second;
            if (p != root && !contains(lnd2par, p)) {
                if (contains(end2par, p)) {
                    parentsOfExits.insert(end2par.at(p));
                    assert(root == nullptr || !parentsOfExits.empty());
                    assert(parentsOfExits.size() < 2);
                } else {
                    assert(p->get_parent() == nullptr);
                }
                root = p;
            }
        }
        if (root == nullptr && firstLaIt == loopEmbeddings.end()) {
            // no loops, the tree is just polytomy for this subproblem
            *treeExpStream << '(';
            OttIdSet ois = shouldHaveBeenPrunedIds;
            bool first = true;
            for (auto pathPtr : childExitForThisTree) {
                if (contains(sc.pruned_subtrees[treeIndex], pathPtr->phylo_child)) {
                    continue;
                }
                auto rids = pathPtr->get_ott_id_set();
                if (rids.size() != 1) {
                    LOG(DEBUG) << "crashing while exporting OTT ID " << scaffold_node.get_ott_id() << " for tree " << treePtr->get_name();
                    db_write_ott_id_set(" rids = ", rids);
                    assert(false);
                }
                OttId ottId = *rids.begin();
                assert(ottId != LONG_MAX);
                ois.insert(ottId);
            }
            totalLeafSet.insert(ois.begin(), ois.end());
            for (auto ottId : ois) {
                if (!first) {
                    *treeExpStream << ',';
                }
                *treeExpStream << "ott" << ottId;
                first = false;
            }
            *treeExpStream << ");\n";
        } else {
            LOG(DEBUG) << "Before adding child exits...";
            for (auto np : lnd2par) {
                LOG(DEBUG) << "   mapping " << np.first << " " << get_designator(*np.first);
                LOG(DEBUG) << "   to par  " << np.second << " " << get_designator(*np.second);
            }
            LOG(DEBUG) << " root = " << reinterpret_cast<long>(root)  << "  parentsOfExits.size() = " << parentsOfExits.size();
            assert(parentsOfExits.empty() || root != nullptr);
            assert(parentsOfExits.size() < 2 );
             // if the node is uncontested, then all exit paths must have the same parent
            U * const deeperNd = (parentsOfExits.empty() ? nullptr : *begin(parentsOfExits));
            std::set<NodeWithSplits *> toDel;
            for (auto n2pEl : lnd2par) {
                if (contains(sc.pruned_subtrees[treeIndex], n2pEl.first)) {
                    toDel.insert(n2pEl.first);
                }
            }
            for (auto td : toDel) {
                lnd2par.erase(td);
            }
            for (auto pp : childExitForThisTree) {
                if (contains(sc.pruned_subtrees[treeIndex], pp->phylo_child)) {
                    continue;
                }
                assert(!contains(lnd2par, pp->phylo_child));
                if (pp->phylo_parent == root) {
                    lnd2par[pp->phylo_child] = root;
                } else {
                    if (deeperNd == pp->phylo_parent) {
                        lnd2par[pp->phylo_child] = root;
                    } else {
                        if (!contains(lnd2par, pp->phylo_parent)) {
                            if (!contains(end2par, pp->phylo_parent)) {
                                LOG(ERROR) << "crashing on treeIndex " << treeIndex;
                                LOG(ERROR) << "deeperNd is set to  " << reinterpret_cast<long>(deeperNd);
                                debugPrintPathPairing(*pp);
                                assert(false);
                            }
                            // necessary if we have node maps to scaffold node and is also
                            //  the child of a polytomy that is compatible with scaffoldnode
                            auto ep = end2par[pp->phylo_parent];
                            lnd2par[pp->phylo_parent] = ep;
                            while (!contains(lnd2par, ep) && contains(end2par, ep)) {
                                auto nep = end2par[ep];
                                lnd2par[ep] = nep;
                                ep = nep;
                            }
                            if (root != nullptr && !contains(lnd2par, root) && contains(end2par, root)) {
                                lnd2par[root] = end2par[root];
                            }
                        }
                        lnd2par[pp->phylo_child] = pp->phylo_parent;
                    }
                }
                auto rids = pp->get_ott_id_set();
                assert(rids.size() == 1);
                OttId ottId = *rids.begin();
                assert(ottId != LONG_MAX);
                totalLeafSet.insert(ottId);
                nd2id[pp->phylo_child] = ottId;
            }
            RootedTreeTopologyNoData toWrite;
            LOG(DEBUG) << "After adding child exits...";
            for (auto np : lnd2par) {
                LOG(DEBUG) << "   mapping " << np.first << " " << get_designator(*np.first);
                LOG(DEBUG) << "   to par  " << np.second << " " << get_designator(*np.second);
            }
            try {
                copy_tree_structure(lnd2par, nd2id, toWrite);
            } catch (const OTCError & ) {
                LOG(ERROR) << "could not construct a valid tree";
                debugPrint(scaffold_node, treeIndex, sn2ne);
                debugPrintNd2Par("lnd2par", lnd2par);
                assert(false);
            }
            sort_children_by_lowest_des_ott_id(toWrite.get_root());
            write_tree_as_newick(*treeExpStream, toWrite);
            *treeExpStream << "\n";
        }
        *provExpStream << treePtr->get_name() << "\n";
    }
    GreedyBandedForest<T, U> gpf{scaffold_node.get_ott_id()};
    gpf.attempt_to_add_grouping(totalLeafSet, EMPTY_SET, 0, 1, &sc);
    for (std::size_t treeIndex = 0 ; treeIndex < sc.num_trees; ++treeIndex) {
        for (auto snc : iter_child(scaffold_node)) {
            assert(snc != nullptr);
        }
    }
    gpf.finish_resolution_of_embedded_clade(scaffold_node, this, &sc);
    if (exportStream == nullptr) {
        provFileStream.close();
        treeFileStream.close();
    }
    //debugPrint(scaffold_node, 7, sn2ne);
    //debug_node_embeddings("leaving export_subproblem_and_resolve", false, sn2ne);
    return retStr;
}


template<typename T, typename U>
void NodeEmbedding<T, U>::resolve_given_uncontested_monophyly(T & scaffold_node,
                                                           SupertreeContextWithSplits & sc) {
    const OttIdSet EMPTY_SET;
    LOG(DEBUG) << "resolve_given_uncontested_monophyly for " << scaffold_node.get_ott_id();
    GreedyBandedForest<T, U> gpf{scaffold_node.get_ott_id()};
    std::set<PathPairing<T, U> *> considered;
    const auto scaffOTTId = scaffold_node.get_ott_id();
    std::string forestDOTfile = "forestForOTT";
    forestDOTfile += std::to_string(scaffOTTId);
    for (std::size_t treeIndex = 0 ; treeIndex < sc.num_trees; ++treeIndex) {
        const auto laIt = loopEmbeddings.find(treeIndex);
        if (laIt == loopEmbeddings.end()) {
            continue;
        }
        LOG(INFO) << "      treeIndex = " << treeIndex;
        const OttIdSet relevantIds = get_relevant_des_ids(sc.scaffold_to_node_embedding, treeIndex);
        PathPairSet & pps = laIt->second;
        // leaf set of this tree for this subtree
        // for repeatability, we'll try to add groupings in reverse order of des_ids sets (deeper first)
        std::map<OttIdSet, PathPairing<T, U> *> mapToProvideOrder;
        for (auto pp : pps) {
            mapToProvideOrder[pp->get_ott_id_set()] = pp;
        }
        OttId bogusGroupIndex = 0; // should get this from the node!
        typedef std::pair<const OttIdSet *, PathPairing<T, U> *>  q_t;
        std::queue<q_t> trivialQ;
        for (auto mpoIt = mapToProvideOrder.rbegin(); mpoIt != mapToProvideOrder.rend(); ++mpoIt) {
            auto ppptr = mpoIt->second;
            if (ppptr->path_is_now_trivial()) {
                db_write_ott_id_set("path_is_now_trivial", mpoIt->first);
                const q_t toQ{&(mpoIt->first), ppptr};
                trivialQ.push(toQ);
            } else {
                const auto & d = mpoIt->first;
                gpf.debug_invariants_check();
                LOG(INFO) << "        bogusGroupIndex = " << bogusGroupIndex << " out of " << mapToProvideOrder.size() << " (some of which may be skipped as trivial)";
                gpf.attempt_to_add_grouping(d, relevantIds, static_cast<int>(treeIndex), bogusGroupIndex, &sc);
                gpf.debug_invariants_check();
                considered.insert(ppptr);
            }
            bogusGroupIndex++;
        }
        while (!trivialQ.empty()) {
            const q_t triv = trivialQ.front();
            const OttIdSet * inc = triv.first;
            const auto ppptr = triv.second;
            gpf.add_leaf(*inc, relevantIds, static_cast<int>(treeIndex), bogusGroupIndex++, &sc);
            gpf.debug_invariants_check();
            considered.insert(ppptr);
            trivialQ.pop();
        }
    }
    // we might have missed some descendants  - any child that is has
    //  "scaffold_node" as its embedded parent, but which is not involved
    //  any loop or exiting edges.
    //  This means that we have no info on the placement of such nodes.
    //      so we'll just attach them here.
    //  First step: get the list of paths for the children.
    std::size_t bogusTreeIndex = 123456; // should get this from the node!
    OttId bogusGroupIndex = 100000; // should get this from the node!
    auto childExitPaths = get_all_child_exit_paths(scaffold_node, sc.scaffold_to_node_embedding);
    for (auto pathPtr : childExitPaths) {
        if (!contains(considered, pathPtr)) {
            gpf.attempt_to_add_grouping(pathPtr->get_ott_id_set(),
                                    EMPTY_SET,
                                    static_cast<int>(bogusTreeIndex),
                                    bogusGroupIndex++,
                                    &sc);
            gpf.debug_invariants_check();
            considered.insert(pathPtr); // @TMP not needed
        }
    }
    for (std::size_t treeIndex = 0 ; treeIndex < sc.num_trees; ++treeIndex) {
        for (auto snc : iter_child(scaffold_node)) {
            assert(snc != nullptr);
        }
    }
    gpf.finish_resolution_of_embedded_clade(scaffold_node, this, &sc);
}

template<typename T, typename U>
std::map<std::size_t, std::set<PathPairing<T, U> *> > copyAllLoopPathPairing(const T *nd, const std::map<const T *, NodeEmbedding<T, U> > & eForNd) {
    const NodeEmbedding<T, U> & ne = eForNd.at(nd);
    return ne.loopEmbeddings;
}

template<typename T, typename U>
void debugPrintPathPairing(const PathPairing<T,U> & clp) {
    std::cerr << "  PathPairing " << reinterpret_cast<long>(&clp) << " scaffold_anc (@" << reinterpret_cast<long>(clp.scaffold_anc) << ") = ott" << clp.scaffold_anc->get_ott_id() << "\n      ";
    //write_newick(std::cerr, clp.scaffold_anc);
    std::cerr << "            scaffold_des (@" << reinterpret_cast<long>(clp.scaffold_des) << ") = ott" << clp.scaffold_des->get_ott_id() << "\n      ";
    //write_newick(std::cerr, clp.scaffold_des);
    std::cerr << "            phylo_parent (@" << reinterpret_cast<long>(clp.phylo_parent) << ") = " << get_designator(*clp.phylo_parent) << "\n      ";
    //write_newick(std::cerr, clp.phylo_parent);
    std::cerr << "            phylo_child (@" << reinterpret_cast<long>(clp.phylo_child) << ") = " << get_designator(*clp.phylo_child) << "\n";
    //write_newick(std::cerr, clp.phylo_child);
    
}

template<typename T, typename U>
void NodeEmbedding<T, U>::debugPrint(T & scaffold_node,
                                     std::size_t treeIndex,
                                     const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const {
    for (auto child : iter_child(scaffold_node)) {
        auto & cne = sn2ne.at(child);
        auto cneIt = cne.edgeBelowEmbeddings.find(treeIndex);
        if (cneIt == cne.edgeBelowEmbeddings.end()) {
            continue;
        }
        std::cerr << " Child ott" << child->get_ott_id() << "\n";
        for (const auto & cneExit : cneIt->second) {
            debugPrintPathPairing(*cneExit);
        }
    }
    auto & sne = sn2ne.at(&scaffold_node);
    std::cerr << " Loops for parent ott" << scaffold_node.get_ott_id() << "\n";
    auto sneIt = sne.loopEmbeddings.find(treeIndex);
    if (sneIt != sne.loopEmbeddings.end()) {
        for (auto & snExit : sneIt->second) {
            debugPrintPathPairing(*snExit);
        }
    }
    std::cerr << " Exits for parent ott" << scaffold_node.get_ott_id() << "\n";
    sneIt = sne.edgeBelowEmbeddings.find(treeIndex);
    if (sneIt != sne.edgeBelowEmbeddings.end()) {
        for (auto & snExit : sneIt->second) {
            debugPrintPathPairing(*snExit);
        }
    }
}

template<typename T, typename U>
void NodeEmbedding<T, U>::prune_suppressed(std::size_t treeIndex, U * phyloPar, U * phylo_child) {
    if (phyloPar->get_out_degree() == 2) {
        auto & pps = loopEmbeddings.at(treeIndex);
        PathPairSet toDel;
        for (auto pp : pps) {
            if (pp->phylo_child == phyloPar) {
                toDel.insert(pp);
            }
        }
        assert(toDel.size() < 2);
        for (auto pp : toDel) {
            pps.erase(pp);
        }
    }
    phyloPar->remove_child(phylo_child);
}

template<typename T, typename U>
void registerAsPruned(std::size_t treeIndex, T * phylo_child, U & sc) {
    auto & ps = sc.pruned_subtrees[treeIndex];
    ps.insert(phylo_child);
    auto p = phylo_child->get_parent();
    for (auto c : iter_child(*p)) {
        if (!contains(ps, c)) {
            return;
        }
    }
    registerAsPruned(treeIndex, p, sc);
}
template<typename T, typename U>
void NodeEmbedding<T, U>::collapse_group(T & scaffold_node, SupertreeContext<T, U> & sc) {
    assert(&scaffold_node == embeddedNode);
    sc.log(COLLAPSE_TAXON, scaffold_node);
    assert(!scaffold_node.is_tip());
    U * p = scaffold_node.get_parent();
    assert(p != nullptr); // can't disagree with the root !
    // remap all nodes in NodePairing to parent
    for (auto nai : nodeEmbeddings) {
        for (auto np : nai.second) {
            np->scaffold_node = p;
        }
    }
    std::map<const T *, NodeEmbedding<T, U> > & sn2ne = sc.scaffold_to_node_embedding;
    NodeEmbedding<T, U> & parEmbedding = const_cast<NodeEmbedding<T, U> &>(sn2ne.at(p));
    LOG(DEBUG) << "TOP of collapse_group";
    //parEmbedding.debugPrint(scaffold_node, 7, sn2ne);
    //const auto beforePL = copyAllLoopPathPairing(p, sn2ne);
    // every loop for this node becomes a loop for its parent
    for (auto lai : loopEmbeddings) {
        const auto & treeIndex = lai.first;
        for (auto lp : lai.second) {
            assert(lp->scaffold_des == &scaffold_node);
            assert(lp->scaffold_anc == &scaffold_node);
            lp->scaffold_des = p;
            lp->scaffold_anc = p;
            parEmbedding.loopEmbeddings[treeIndex].insert(lp);
        }
    }
    std::set<std::size_t> indsOfTreesWithNewLoops;
    std::set<std::size_t> indsOfTreesMappedToInternal;
    // every exit edge for this node becomes a loop for its parent if it is not trivial
    for (auto ebai : edgeBelowEmbeddings) {
        const auto & treeIndex = ebai.first;
        std::set<PathPairPtr> pathsAblated;
        for (auto lp : ebai.second) {
            if (lp->phylo_child->is_tip()
                && lp->scaffold_des == &scaffold_node) {
                // this only happens if a terminal was mapped to this higher level taxon
                // we don't know how to interpret this label any more, so we should drop that 
                // leaf. Unless specifically requested not to by the user.
                // The taxa will be included by other relationships (the taxonomy as
                // a last resort), so we don't need to worry about losing leaves by skipping this...
                if (sc.prune_tips_mapped_to_contested_taxa) {
                    LOG(INFO) << "Ablating path leading to ott" << lp->phylo_child->get_ott_id();;
                    LOG(DEBUG) << "IGNORING scaff = " << scaffold_node.get_ott_id() << " == phylo " << lp->phylo_child->get_ott_id();
                    assert(scaffold_node.get_ott_id() == lp->phylo_child->get_ott_id());
                    sc.log(IGNORE_TIP_MAPPED_TO_NONMONOPHYLETIC_TAXON, *lp->phylo_child);
                    OttIdSet innerOTTId;
                    innerOTTId.insert(lp->phylo_child->get_ott_id());
                    OttIdSet n = scaffold_node.get_data().des_ids; // expand the internal name to it taxonomic content
                    n.erase(lp->phylo_child->get_ott_id());
                    lp->update_des_ids_for_self_and_anc(innerOTTId, n, sn2ne);
                    indsOfTreesMappedToInternal.insert(treeIndex);
                    pathsAblated.insert(lp);
                    registerAsPruned(treeIndex, lp->phylo_child, sc);
                } else {
                    phyloNd2ParForUnembeddedTrees[treeIndex][lp->phylo_child] = lp->phylo_parent;
                }
            } else if (lp->scaffold_anc == p) {
                //if (lp->scaffold_des == &scaffold_node) {
                if ((!lp->phylo_child->is_tip()) && lp->scaffold_des == &scaffold_node) {
                    lp->scaffold_des = p;
                    parEmbedding.loopEmbeddings[treeIndex].insert(lp);
                    indsOfTreesWithNewLoops.insert(treeIndex);
                }
                /*} else {
                    if (lp->scaffold_des->get_parent() != &scaffold_node) {
                        LOG(ERROR) << " Anbandonding a path that ends at ott" << lp->scaffold_des->get_ott_id();
                        LOG(ERROR) << "                    and starts at ott" << lp->scaffold_anc->get_ott_id() << " when collapsing ott" << scaffold_node.get_ott_id();
                        LOG(ERROR) << " lp->scaffold_des->get_parent() ott" << lp->scaffold_des->get_parent()->get_ott_id();
                        assert(false);
                    }
                }*/
            } else {
                // if the anc isn't the parent, then it must pass through scaffold_node's par
                assert(contains(parEmbedding.edgeBelowEmbeddings[treeIndex], lp));
                if (lp->scaffold_des == &scaffold_node && !lp->phylo_child->is_tip()) {
                    lp->scaffold_des = p;
                }
            }
        }
        if (!pathsAblated.empty()) {
            for (auto deadPathPtr : pathsAblated) {
                debugPrintPathPairing(*deadPathPtr);
                /*T * ablatedScaffAnc = deadPathPtr->scaffold_anc;
                U * ablatedPhyloChild = deadPathPtr->phylo_child;
                U * ablatedPhyloPar = deadPathPtr->phylo_parent;
                U * curr = &scaffold_node;
                while (curr != deadPathPtr->scaffold_anc) {
                    sn2ne.at(curr).remove_ref_to_exit_path(treeIndex, deadPathPtr);
                    curr = curr->get_parent();
                }
                //sn2ne.at(ablatedScaffAnc).prune_suppressed(treeIndex, ablatedPhyloPar, ablatedPhyloChild);
                */
            }
        }
    }
    for (auto child : iter_child(scaffold_node)) {
        auto cit = sn2ne.find(child);
        if (cit == sn2ne.end()) {
            continue;
        }
        NodeEmbedding<T, U>& childEmbedding = cit->second;
        for (auto ceabi : childEmbedding.edgeBelowEmbeddings) {
            for (auto clp : ceabi.second) {
                if (clp->scaffold_anc == &scaffold_node) {
                    clp->scaffold_anc = p;
                }
            }
        }
    }
    prune_collapsed_node(scaffold_node, sc);
    /*
    const auto afterPL = copyAllLoopPathPairing(p, sn2ne);

    for (auto bpl : beforePL) {
        const auto & afterVal = afterPL.at(bpl.first);
        const auto tv = parEmbedding.get_all_incoming_path_pairs(sn2ne, bpl.first);
        assert(is_subset(bpl.second, afterVal));
        for (auto pp : bpl.second) {
            const PathPairing<T, U> * cpp = pp;
            LOG(DEBUG) << "Checking tree" << bpl.first << " scaff" << p->get_ott_id() << " " << cpp->scaffold_anc->get_ott_id() << " -> " <<  cpp->scaffold_des->get_ott_id();
            assert(vcontains(tv, cpp));
        }
    }
    */
}

template<typename T, typename U>
void NodeEmbedding<T, U>::prune_collapsed_node(T & scaffold_node, SupertreeContextWithSplits & sc) {
    check_all_node_pointers_iter(scaffold_node);
     LOG(DEBUG) << "collapsed paths from ott" << scaffold_node.get_ott_id() << ", adding child to parent";
    // NOTE: it is important that we add the children of scaffold_node the left of its location
    //  in the tree so that the postorder traversal will not iterate over them.
    assert(scaffold_node.has_children());
    auto p = scaffold_node.get_parent();
    assert(p);
    if (!phyloNd2ParForUnembeddedTrees.empty()) {
        auto & scaffoldPar = sc.scaffold_to_node_embedding.at(p);
        // We have a hacky data structure to pass along to our parent if
        //  we were asked to retain terminal nodes that are mapped to this
        //  these tips will no longer be embedded in the scaffold... Ugh.
        for (auto pmfu : phyloNd2ParForUnembeddedTrees) {
            if (pmfu.second.empty()) {
                continue;
            }
            const auto treeIndex = pmfu.first;
            auto & ppmfuForTree =  scaffoldPar.phyloNd2ParForUnembeddedTrees[treeIndex];
            for (auto dpmfuForTree : pmfu.second) {
                ppmfuForTree[dpmfuForTree.first] = dpmfuForTree.second;
            }
        }
    }
    while(scaffold_node.has_children())
    {
        auto n = scaffold_node.get_first_child();
        n->detach_this_node();
        scaffold_node.add_sib_on_left(n);
    }
    scaffold_node.detach_this_node();
    sc.scaffold_tree.mark_as_detached(&scaffold_node);
    sc.detached_scaffold_nodes.insert(&scaffold_node);
}

template<typename T, typename U>
OttIdSet NodeEmbedding<T, U>::get_relevant_des_ids(const std::map<const T *, NodeEmbedding<T, U> > & eForNd,
                                                std::size_t treeIndex) {
    /* find MRCA of the phylo nodes */
    auto ippV = get_all_incoming_path_pairs(eForNd, treeIndex);
    OttIdSet relevantIds;
    for (auto pIt : ippV) {
        const OttIdSet otherRelevantIds = get_relevant_des_ids_from_path(*pIt);
        relevantIds.insert(otherRelevantIds.begin(), otherRelevantIds.end());
    }
    return relevantIds;
}

template<typename T, typename U>
bool NodeEmbedding<T, U>::report_if_contested(std::ostream & out,
                       const T * nd,
                       const std::vector<TreeMappedWithSplits *> & treePtrByIndex,
                       const std::vector<NodeWithSplits *> & aliasedBy,
                       bool verbose) const {
    if (is_contested()) {
        auto c = get_contesting_tree_indices();
        for (auto cti : c) {
            auto ctree = treePtrByIndex.at(cti);
            const OttIdSet ls = get_ott_id_set_for_leaves(*ctree);
            const auto & edges = get_edges_exiting(cti);
            const std::string prefix = get_contested_preamble(*nd, *ctree);
            if (verbose) {
                report_on_conflicting(out, prefix, nd, edges, ls);
            } else {
                out << prefix << '\n';
            }
            for (auto na : aliasedBy) {
                const std::string p2 = get_contested_preamble(*na, *ctree);
                if (verbose) {
                    report_on_conflicting(out, p2, na, edges, ls);
                } else {
                    out << prefix << '\n';
                }
            }
        }
        return true;
    }
    return false;
}

template<typename T, typename U>
void report_on_conflicting(std::ostream & out,
                        const std::string & prefix,
                        const T * scaffold,
                        const std::set<PathPairing<T, U> *> & exitPaths,
                        const OttIdSet & phyloLeafSet) {
    if (exitPaths.size() < 2) {
        assert(false);
        throw OTCError("asserts are disabled, but one is not true");
    }
    const auto scaffold_des = set_intersection_as_set(scaffold->get_data().des_ids, phyloLeafSet);
    auto epIt = begin(exitPaths);
    const PathPairing<T, U> * ep = *epIt;
    const U * phyloPar = ep->phylo_parent;
    const U * deepestPhylo = nullptr;
    std::map<OttIdSet, const U *> desIdSet2NdConflicting;
    if (is_proper_subset(scaffold_des, phyloPar->get_data().des_ids)) {
        deepestPhylo = phyloPar;
    } else {
        desIdSet2NdConflicting[phyloPar->get_data().des_ids] = phyloPar;
        for (auto anc : iter_anc_const(*phyloPar)) {
            if (is_proper_subset(scaffold_des, anc->get_data().des_ids)) {
                deepestPhylo = anc;
                break;
            }
            desIdSet2NdConflicting[anc->get_data().des_ids] = anc;
        }
        assert(deepestPhylo != nullptr);
    }
    for (++epIt; epIt != end(exitPaths); ++epIt) {
        const U * phyloNd  = (*epIt)->phylo_child;
        assert(phyloNd != nullptr);
        for (auto anc : iter_anc_const(*phyloNd)) {
            if (anc == deepestPhylo) {
                break;
            }
            desIdSet2NdConflicting[anc->get_data().des_ids] = anc;
        }
    }
    if (desIdSet2NdConflicting.empty()) {
        out << "VERY ODD " << prefix << ", desIdSet2NdConflicting but desIdSet2NdConflicting is empty()!\n";
        //assert(false); // not reachable if we are calling is_contested first as a test.
        throw OTCError("Expecting is_contested to have guarded call to report_on_conflicting");
    }
    for (const auto & mIt : desIdSet2NdConflicting) {
        const auto & di = mIt.first;
        auto nd = mIt.second;
        const OttIdSet e = set_difference_as_set(di, scaffold_des);
        const OttIdSet m = set_difference_as_set(scaffold_des, di);
        out << prefix;
        emit_conflict_details(out, *nd, e, m);
    }
}

template class NodePairing<NodeWithSplits, NodeWithSplits>; // force explicit instantiaion of this template.
template class PathPairing<NodeWithSplits, NodeWithSplits>; // force explicit instantiaion of this template.
template class NodeEmbedding<NodeWithSplits, NodeWithSplits>; // force explicit instantiaion of this template.

}// namespace
