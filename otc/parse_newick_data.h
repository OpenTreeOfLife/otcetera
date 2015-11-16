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

void newickParseNodeInfo(RootedTree<RTNodeNoData, RTreeNoData> & ,
                         RootedTreeNode<RTNodeNoData> &,
                         const NewickTokenizer::Token * labelToken,
                         const NewickTokenizer::Token * colonToken, // used for comment
                         const NewickTokenizer::Token * branchLenToken,
                         const ParsingRules & parsingRules);

////////////////////////////////////////////////////////////////////////////////
// no ops for no data
inline void newickCloseNodeHook(RootedTree<RTNodeNoData, RTreeNoData> & ,
                                RootedTreeNode<RTNodeNoData> & ,
                                const NewickTokenizer::Token & ,
                                const ParsingRules & ) {
}

inline void newickParseNodeInfo(RootedTree<RTNodeNoData, RTreeNoData> & ,
                                RootedTreeNode<RTNodeNoData> & node,
                                const NewickTokenizer::Token * labelToken,
                                const NewickTokenizer::Token * , // used for comment
                                const NewickTokenizer::Token * ,
                                const ParsingRules &parsingRules) {
    if (labelToken) {
        node.setName(labelToken->content());
        if (not parsingRules.setOttIds)
            return;
        if ((!parsingRules.setOttIdForInternals) && node.isInternal()) {
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
            node.setOttId(ottID);
        }
    }
}

template<typename N, typename T>
inline void postParseHook(RootedTree<N,T> & , const ParsingRules & ) {
}
////////////////////////////////////////////////////////////////////////////////
// tree-level mapping of ottID to Node data
inline void newickCloseNodeHook(RootedTree<RTNodeNoData, RTreeOttIDMapping<RTNodeNoData> > & ,
                                RootedTreeNode<RTNodeNoData> & ,
                                const NewickTokenizer::Token & ,
                                const ParsingRules & ) {
}

template <typename T, typename U>
void setOttIdAndAddToMap(RootedTree<T, U> & tree,
                         RootedTreeNode<T> & node,
                         const NewickTokenizer::Token * labelToken,
                         const ParsingRules &);

template <typename T, typename U>
inline void setOttIdAndAddToMap(RootedTree<T, U> & tree,
                         RootedTreeNode<T> & node,
                         const NewickTokenizer::Token * labelToken,
                         const ParsingRules & parsingRules) {
    if (labelToken) {
        node.setName(labelToken->content());
        if (not parsingRules.setOttIds) return;
        if (not parsingRules.setOttIdForInternals and node.isInternal()) return;
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
            node.setOttId(ottID);
            U & treeData = tree.getData();
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

inline void newickParseNodeInfo(RootedTree<RTNodeNoData, RTreeOttIDMapping<RTNodeNoData> > & tree,
                                RootedTreeNode<RTNodeNoData> & node,
                                const NewickTokenizer::Token * labelToken,
                                const NewickTokenizer::Token * , // used for comment
                                const NewickTokenizer::Token * ,
                                const ParsingRules & parsingRules) {
    setOttIdAndAddToMap(tree, node, labelToken, parsingRules);
}

////////////////////////////////////////////////////////////////////////////////
// tree-level mapping of ottID to Node data
inline void newickCloseNodeHook(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & ,
                                RootedTreeNode<RTSplits> & node,
                                const NewickTokenizer::Token & token,
                                const ParsingRules & parsingRules) {
    if (not parsingRules.setOttIds) return;
    if (node.isTip()
        && !node.hasOttId()
        && (!parsingRules.pruneUnrecognizedInputTips)) {
        throw OTCParsingContentError("Expecting each tip to have an ID.", token.getStartPos());
    }
}

inline void newickParseNodeInfo(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & tree,
                                RootedTreeNode<RTSplits> & node,
                                const NewickTokenizer::Token * labelToken,
                                const NewickTokenizer::Token * , // used for comment
                                const NewickTokenizer::Token * ,
                                const ParsingRules & parsingRules) {
    setOttIdAndAddToMap(tree, node, labelToken, parsingRules);
}

template<>
inline void postParseHook(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & tree, const ParsingRules & parsingRules) {
    if (not parsingRules.setOttIds) return;
    if (parsingRules.includeInternalNodesInDesIdSets) {
        fillDesIdSetsIncludingInternals(tree);
    } else {
        fillDesIdSets(tree);
    }
    for (auto p : tree.getData().ottIdToNode) {
        assert(p.first == p.second->getOttId());
    }
}

} // namespace otc

#endif
