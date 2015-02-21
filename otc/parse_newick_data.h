#ifndef OTCETERA_PARSE_NEWICK_DATA_H
#define OTCETERA_PARSE_NEWICK_DATA_H

#include "otc/otc_base_includes.h"
#include "otc/tree.h"
#include "otc/tree_data.h"
#include "otc/newick_tokenizer.h"
#include "otc/tree_iter.h"

namespace otc {

void newickParseNodeInfo(RootedTree<RTNodeNoData, RTreeNoData> & ,
						 RootedTreeNode<RTNodeNoData> &,
						 const NewickTokenizer::Token * labelToken,
						 const NewickTokenizer::Token * colonToken, // used for comment
						 const NewickTokenizer::Token * branchLenToken);
void newickCloseNodeHook(RootedTree<RTNodeNoData, RTreeNoData> & ,
						 RootedTreeNode<RTNodeNoData> &,
						 const NewickTokenizer::Token & closeToken);
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
		node.SetName(labelToken->content());
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
			node.SetOTUId(ottID);
			RTreeOttIDMapping<RTNodeNoData> & treeData = tree.GetData();
			treeData.ottIdToNode[ottID] = &node;
		}
	}
}
inline void postParseHook(RootedTree<RTNodeNoData, RTreeOttIDMapping<RTNodeNoData> > & ) {
}
////////////////////////////////////////////////////////////////////////////////
// tree-level mapping of ottID to Node data
inline void newickCloseNodeHook(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & ,
								RootedTreeNode<RTSplits> & ,
								const NewickTokenizer::Token & ) {
}

inline void newickParseNodeInfo(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & tree,
								RootedTreeNode<RTSplits> & node,
								const NewickTokenizer::Token * labelToken,
								const NewickTokenizer::Token * , // used for comment
								const NewickTokenizer::Token * ) {
	if (labelToken) {
		long ottID = ottIDFromName(labelToken->content());
		if (ottID >= 0) {
			node.SetOTUId(ottID);
			RTreeOttIDMapping<RTSplits> & treeData = tree.GetData();
			treeData.ottIdToNode[ottID] = &node;
		}
	}
}

inline void postParseHook(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & tree) {
	for (auto node : PostorderInternalNode<RTSplits, RTreeOttIDMapping<RTSplits> >(tree)) {
		std::set<long> & mrca = node->GetData().mrca;
		if (node->IsTip()) {
			mrca.insert(node->GetOTUId());
		} else {
			for (auto child : ChildIterator<RTSplits, RTreeOttIDMapping<RTSplits> >(*node)) {
				std::set<long> & cmrca = child->GetData().mrca;
				mrca.insert(cmrca.begin(), cmrca.end());
			}
		}
	}
}

} // namespace otc

#endif
