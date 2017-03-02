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
inline void newick_close_node_hook(RootedTree<N,T> & ,
                                RootedTreeNode<N> & ,
                                const NewickTokenizer::Token & ,
                                const ParsingRules & )
{ }

template <typename N, typename T>
inline void newick_parse_node_info(RootedTree<N, T> & ,
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
        long ottID = ott_id_from_name(labelToken->content());
        if (ottID >= 0) {
            if (parsingRules.id_remapping != nullptr) {
                auto rIt = parsingRules.id_remapping->find(ottID);
                if (rIt != parsingRules.id_remapping->end()) {
                    LOG(DEBUG) << "id_remapping from OTT" << ottID << " to OTT" << rIt->second;
                    ottID = rIt->second;
                }
            }
            node.set_ott_id(ottID);
        }
    }
}

template<typename N, typename T>
inline void post_parse_hook(RootedTree<N,T> & , const ParsingRules & ) {
}
////////////////////////////////////////////////////////////////////////////////
// tree-level mapping of ottID to Node data

template <typename T, typename U>
void set_ott_id_and_add_to_map(RootedTree<T, U> & tree,
                         RootedTreeNode<T> & node,
                         const NewickTokenizer::Token * labelToken,
                         const ParsingRules &);

template <typename T, typename U>
inline void set_ott_id_and_add_to_map(RootedTree<T, U> & tree,
                         RootedTreeNode<T> & node,
                         const NewickTokenizer::Token * labelToken,
                         const ParsingRules & parsingRules) {
    if (labelToken) {
        node.setName(labelToken->content());
        if (not parsingRules.set_ott_ids) return;
        if (not parsingRules.set_ott_idForInternals and node.isInternal()) return;
        long ottID = ott_id_from_name(labelToken->content());
        if (ottID >= 0) {
            if (parsingRules.ott_id_validator != nullptr) {
                if (!contains(*parsingRules.ott_id_validator, ottID)) {
                    if (!parsingRules.prune_unrecognized_input_tips) {
                        std::string m = "Unrecognized OTT Id ";
                        m += std::to_string(ottID);
                        throw OTCParsingError(m.c_str(), labelToken->content(), labelToken->get_start_pos());
                    } else {
                        return;
                    }
                }
            }
            if (parsingRules.id_remapping != nullptr) {
                auto rIt = parsingRules.id_remapping->find(ottID);
                if (rIt != parsingRules.id_remapping->end()) {
                    LOG(DEBUG) << "id_remapping from OTT" << ottID << " to OTT" << rIt->second;
                    ottID = rIt->second;
                }
            }
            node.set_ott_id(ottID);
            U & treeData = tree.get_data();
            if (contains(treeData.ottIdToNode, ottID)) {
                throw OTCParsingError("Expecting an OTT Id to only occur one time in a tree.",
                                      labelToken->content(),
                                      labelToken->get_start_pos());
            }
            treeData.ottIdToNode[ottID] = &node;
        } else if (parsingRules.require_ott_ids and not parsingRules.prune_unrecognized_input_tips) {
            throw OTCParsingError("Expecting a name for a taxon to end with an ott##### where the numbers are the OTT Id.",
                                  labelToken->content(),
                                  labelToken->get_start_pos());
        }
    }
}

template<typename N>
inline void newick_parse_node_info(RootedTree<N, RTreeOttIDMapping<N> > & tree,
                                RootedTreeNode<N> & node,
                                const NewickTokenizer::Token * labelToken,
                                const NewickTokenizer::Token * , // used for comment
                                const NewickTokenizer::Token * ,
                                const ParsingRules & parsingRules) {
    set_ott_id_and_add_to_map(tree, node, labelToken, parsingRules);
}

////////////////////////////////////////////////////////////////////////////////
// tree-level mapping of ottID to Node data
template<>
inline void newick_close_node_hook(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & ,
                                RootedTreeNode<RTSplits> & node,
                                const NewickTokenizer::Token & token,
                                const ParsingRules & parsingRules) {
    if (not parsingRules.set_ott_ids) return;
    if (node.isTip()
        && !node.has_ott_id()
        && (!parsingRules.prune_unrecognized_input_tips)) {
        throw OTCParsingContentError("Expecting each tip to have an ID.", token.get_start_pos());
    }
}

template<>
inline void post_parse_hook(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & tree, const ParsingRules & parsingRules) {
    if (not parsingRules.set_ott_ids) return;
    if (parsingRules.include_internal_nodes_in_des_id_sets) {
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
