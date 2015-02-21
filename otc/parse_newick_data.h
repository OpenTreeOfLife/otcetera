#ifndef OTCETERA_PARSE_NEWICK_DATA_H
#define OTCETERA_PARSE_NEWICK_DATA_H

#include "otc/otc_base_includes.h"
#include "otc/tree.h"
#include "otc/tree_data.h"
#include "otc/newick_tokenizer.h"

namespace otc {




void newickParseNodeInfo(RootedTree<RTNodeNoData, RTreeNoData> & ,
						 RootedTreeNode<RTNodeNoData> &,
						 const NewickTokenizer::Token * labelToken,
						 const NewickTokenizer::Token * colonToken, // used for comment
						 const NewickTokenizer::Token * branchLenToken);
void newickCloseNodeHook(RootedTree<RTNodeNoData, RTreeNoData> & ,
						 RootedTreeNode<RTNodeNoData> &,
						 const NewickTokenizer::Token & closeToken);

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
		node.SetName(labelToken->content());
	}
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
			node.SetOTUId(ottID);
			RTreeOttIDMapping<RTNodeNoData> & treeData = tree.GetData();
			treeData.ottIdToNode[ottID] = &node;
		}
	}
}

} // namespace otc

#endif
