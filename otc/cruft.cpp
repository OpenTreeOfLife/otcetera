

// addIngroupAtNode adds leaves for all "new" ottIds to:
//      the root if they are in the excGroup only, OR
//      ing if they are in the incGroup.
// The caller has guaranteed that:
//   outg is a part of the only tree in the forest that intersects with leafSet, AND
//   outg is the MRCA of all of the taxa in leafset in this tree.
//   ing is the MRCA of the incGroup in this tree
//
template<typename T, typename U>
void GreedyPhylogeneticForest<T,U>::addIngroupAtNode(NodeWithSplits *, //delete root param?
                                                     NodeWithSplits *ing,
                                                     NodeWithSplits *outg,
                                                     const OttIdSet & incGroup,
                                                     const OttIdSet & leafSet) {
    assert(outg != nullptr);
    if (ing == nullptr) {
        // there was no intersection with the incGroup...
        // So: create a new node at the MRCA of the excGroup, and add the incGroup to that node.
        assert(outg != nullptr);
        ing = nodeSrc.createChild(outg);
        ing->getData().desIds = incGroup;
        for (auto io : incGroup) {
            auto nn = addChildForOttId(*ing, io, nodeSrc);
            registerTreeForNode(nn, )

        }
    } else {
        const auto newIngLeaves = set_difference_as_set(incGroup, ottIdSet);
        if (!newIngLeaves.empty()) {
            ing->getData().desIds.insert(begin(newIngLeaves), end(newIngLeaves));
            for (auto io : newIngLeaves) {
                addChildForOttId(*ing, io, nodeSrc);
            }
            for (auto ianc : iter_anc(*ing)) {
                if (ianc == outg) {
                    break; // we'll update outg and below in the code below
                }
                ianc->getData().desIds.insert(begin(newIngLeaves), end(newIngLeaves));
            }
        }
    }
    // attach any new excGroup leaves as children of outg, creating/expanding a polytomy.
    const auto newLeaves = set_difference_as_set(leafSet, ottIdSet);
    ottIdSet.insert(begin(newLeaves), end(newLeaves));
    const auto newOutLeaves = set_difference_as_set(newLeaves, incGroup);
    for (auto oo : newOutLeaves) {
        addChildForOttId(*outg, oo, nodeSrc);
    }
    // add new ottids to desIDs in outg and its ancestors
    auto f = outg->getFirstChild();
    assert(f != nullptr);
    for (auto oanc : iter_anc(*f)) {
        oanc->getData().desIds.insert(begin(newLeaves), end(newLeaves));
    }
}

template<typename T, typename U>
void GreedyPhylogeneticForest<T,U>::graftTreesTogether(NodeWithSplits *, //rr,
                                                       NodeWithSplits *, //ri,
                                                       NodeWithSplits *, //delr,
                                                       NodeWithSplits *, //deli,
                                                       NodeWithSplits *, //delo,
                                                       const OttIdSet & , //incGroup,
                                                       const OttIdSet & ){  //leafSet
    NOT_IMPLEMENTED; // refactoring to banded
}

