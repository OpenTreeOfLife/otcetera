#ifndef OTCETERA_PARSE_NEWICK_DATA_H
#define OTCETERA_PARSE_NEWICK_DATA_H

#include "otc/otc_base_includes.h"
#include "otc/util.h"
#include "otc/tree.h"
#include "otc/tree_data.h"
#include "otc/newick_tokenizer.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"

namespace otc {

////////////////////////////////////////////////////////////////////////////////
// no ops for no data
template <typename N, typename T>
inline void newickCloseNodeHook(RootedTree<N,T> & ,
                                RootedTreeNode<N> & ,
                                const NewickTokenizer::Token & ,
                                const ParsingRules & )
{ }

template <typename N, typename T>
inline void newickParseNodeInfo(RootedTree<N, T> & ,
                                RootedTreeNode<N> & node,
                                const NewickTokenizer::Token * labelToken,
                                const NewickTokenizer::Token * , // used for comment
                                const NewickTokenizer::Token * ,
                                const ParsingRules &parsingRules) {
    if (labelToken) {
        node.setName(labelToken->content());
        if (not parsingRules.set_ott_ids)
            return;
        if ((!parsingRules.set_ott_idForInternals) && node.isInternal()) {
            return;
        }
        long ottID = ottIDFromName(labelToken->content());
        if (ottID >= 0) {
            if (parsingRules.idRemapping != nullptr) {
                auto rIt = parsingRules.idRemapping->find(ottID);
                if (rIt != parsingRules.idRemapping->end()) {
                    LOG(DEBUG) << "idRemapping from OTT" << ottID << " to OTT" << rIt->second;
                    ottID = rIt->second;
                }
            }
            node.set_ott_id(ottID);
        }
    }
}

template<typename N, typename T>
inline void postParseHook(RootedTree<N,T> & , const ParsingRules & ) {
}
////////////////////////////////////////////////////////////////////////////////
// tree-level mapping of ottID to Node data

template <typename T, typename U>
void set_ott_idAndAddToMap(RootedTree<T, U> & tree,
                         RootedTreeNode<T> & node,
                         const NewickTokenizer::Token * labelToken,
                         const ParsingRules &);

template <typename T, typename U>
inline void set_ott_idAndAddToMap(RootedTree<T, U> & tree,
                         RootedTreeNode<T> & node,
                         const NewickTokenizer::Token * labelToken,
                         const ParsingRules & parsingRules) {
    if (labelToken) {
        node.setName(labelToken->content());
        if (not parsingRules.set_ott_ids) return;
        if (not parsingRules.set_ott_idForInternals and node.isInternal()) return;
        long ottID = ottIDFromName(labelToken->content());
        if (ottID >= 0) {
            if (parsingRules.ottIdValidator != nullptr) {
                if (!contains(*parsingRules.ottIdValidator, ottID)) {
                    if (!parsingRules.pruneUnrecognizedInputTips) {
                        std::string m = "Unrecognized OTT Id ";
                        m += std::to_string(ottID);
                        throw OTCParsingError(m.c_str(), labelToken->content(), labelToken->getStartPos());
                    } else {
                        return;
                    }
                }
            }
            if (parsingRules.idRemapping != nullptr) {
                auto rIt = parsingRules.idRemapping->find(ottID);
                if (rIt != parsingRules.idRemapping->end()) {
                    LOG(DEBUG) << "idRemapping from OTT" << ottID << " to OTT" << rIt->second;
                    ottID = rIt->second;
                }
            }
            node.set_ott_id(ottID);
            U & treeData = tree.get_data();
            if (contains(treeData.ottIdToNode, ottID)) {
                throw OTCParsingError("Expecting an OTT Id to only occur one time in a tree.",
                                      labelToken->content(),
                                      labelToken->getStartPos());
            }
            treeData.ottIdToNode[ottID] = &node;
        } else if (parsingRules.requireOttIds and not parsingRules.pruneUnrecognizedInputTips) {
            throw OTCParsingError("Expecting a name for a taxon to end with an ott##### where the numbers are the OTT Id.",
                                  labelToken->content(),
                                  labelToken->getStartPos());
        }
    }
}

template<typename N>
inline void newickParseNodeInfo(RootedTree<N, RTreeOttIDMapping<N> > & tree,
                                RootedTreeNode<N> & node,
                                const NewickTokenizer::Token * labelToken,
                                const NewickTokenizer::Token * , // used for comment
                                const NewickTokenizer::Token * ,
                                const ParsingRules & parsingRules) {
    set_ott_idAndAddToMap(tree, node, labelToken, parsingRules);
}

////////////////////////////////////////////////////////////////////////////////
// tree-level mapping of ottID to Node data
template<>
inline void newickCloseNodeHook(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & ,
                                RootedTreeNode<RTSplits> & node,
                                const NewickTokenizer::Token & token,
                                const ParsingRules & parsingRules) {
    if (not parsingRules.set_ott_ids) return;
    if (node.isTip()
        && !node.has_ott_id()
        && (!parsingRules.pruneUnrecognizedInputTips)) {
        throw OTCParsingContentError("Expecting each tip to have an ID.", token.getStartPos());
    }
}

template<>
inline void postParseHook(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & tree, const ParsingRules & parsingRules) {
    if (not parsingRules.set_ott_ids) return;
    if (parsingRules.includeInternalNodesInDesIdSets) {
        fillDesIdSetsIncludingInternals(tree);
    } else {
        fillDesIdSets(tree);
    }
    for (auto p : tree.get_data().ottIdToNode) {
        assert(p.first == p.second->get_ott_id());
    }
}

} // namespace otc

#endif
