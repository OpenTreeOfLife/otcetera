#ifndef OTCETERA_PARSE_NEWICK_DATA_H
#define OTCETERA_PARSE_NEWICK_DATA_H

#include "otc/otc_base_includes.h"
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
						 const NewickTokenizer::Token * branchLenToken);
void postParseHook(RootedTree<RTNodeNoData, RTreeNoData> & );

////////////////////////////////////////////////////////////////////////////////
// no ops for no data
inline void newickCloseNodeHook(RootedTree<RTNodeNoData, RTreeNoData> & ,
								RootedTreeNode<RTNodeNoData> & ,
								const NewickTokenizer::Token & ) {
}

inline void newickParseNodeInfo(RootedTree<RTNodeNoData, RTreeNoData> & ,
								RootedTreeNode<RTNodeNoData> & node,
								const NewickTokenizer::Token * labelToken,
								const NewickTokenizer::Token * , // used for comment
								const NewickTokenizer::Token * ) {
	if (labelToken) {
		node.setName(labelToken->content());
	}
}
inline void postParseHook(RootedTree<RTNodeNoData, RTreeNoData> & ) {
}
////////////////////////////////////////////////////////////////////////////////
// tree-level mapping of ottID to Node data
inline void newickCloseNodeHook(RootedTree<RTNodeNoData, RTreeOttIDMapping<RTNodeNoData> > & ,
								RootedTreeNode<RTNodeNoData> & ,
								const NewickTokenizer::Token & ) {
}

inline void newickParseNodeInfo(RootedTree<RTNodeNoData, RTreeOttIDMapping<RTNodeNoData> > & tree,
								RootedTreeNode<RTNodeNoData> & node,
								const NewickTokenizer::Token * labelToken,
								const NewickTokenizer::Token * , // used for comment
								const NewickTokenizer::Token * ) {
	if (labelToken) {
		long ottID = ottIDFromName(labelToken->content());
		if (ottID >= 0) {
			node.setOttId(ottID);
			RTreeOttIDMapping<RTNodeNoData> & treeData = tree.getData();
			treeData.ottIdToNode[ottID] = &node;
		}
	}
}
inline void postParseHook(RootedTree<RTNodeNoData, RTreeOttIDMapping<RTNodeNoData> > & ) {
}
////////////////////////////////////////////////////////////////////////////////
// tree-level mapping of ottID to Node data
inline void newickCloseNodeHook(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & ,
								RootedTreeNode<RTSplits> & node,
								const NewickTokenizer::Token & token) {
	if (!node.hasOttId()) {
		throw OTCParsingContentError("Expecting each tip to have an ID.", token.getStartPos());
	}
}

inline void newickParseNodeInfo(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & tree,
								RootedTreeNode<RTSplits> & node,
								const NewickTokenizer::Token * labelToken,
								const NewickTokenizer::Token * , // used for comment
								const NewickTokenizer::Token * ) {
	if (labelToken) {
		long ottID = ottIDFromName(labelToken->content());
		if (ottID >= 0) {
			node.setOttId(ottID);
			RTreeOttIDMapping<RTSplits> & treeData = tree.getData();
			treeData.ottIdToNode[ottID] = &node;
		}
	}
}

inline void postParseHook(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & tree) {
	fillDesIdSets<RTSplits, RTreeOttIDMapping<RTSplits> >(tree);
}

} // namespace otc

#endif
