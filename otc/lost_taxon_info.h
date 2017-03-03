#ifndef OTCETERA_LOST_TAXON_INFO_H
#define OTCETERA_LOST_TAXON_INFO_H
#include <map>
#include <iostream>
#include <set>
#include <stack>
#include <vector>
#include "json.hpp"
#include "otc/otc_base_includes.h"
#include "otc/util.h"

namespace otc {
// Stores pointers for a taxon X which not in the solution.
//   points to the MRCA, and all of the attachment points for 
//   subtrees that are part of this taxon.
// An attachment point is:
//      1. an ancestor of at least one member of the taxon X,
//      2. an ancestor of at least one taxon that is not contained in the taxon X.
//      3. has at least one child (an attached node) that consists entirely
//          of member of taxon X
template<typename N>
class LostTaxonDetails {
    using node_vec = std::set<const N *>;
    using attachment_to_node_vec = std::map<const N *, node_vec>;
    public:
    LostTaxonDetails(long oid=0, const N *mrca_nd=nullptr)
        :ott_id(oid),
        mrca(mrca_nd) {
    }
    void addAttachmentPoint(const N * attach, const N * child) {
        attach_node_to_attached_vec[attach].insert(child);
    }
    const N * get_mrca() const {
        return mrca;
    }
    const attachment_to_node_vec & get_attach_node_to_attached_vec() const {
        return attach_node_to_attached_vec;
    }
    attachment_to_node_vec & get_attach_node_to_attached_vec() {
        return attach_node_to_attached_vec;
    }
    const node_vec & get_attached_vec(const N * n) const {
        return attach_node_to_attached_vec.at(n);
    }
    node_vec & get_attached_vec(const N * n) {
        return attach_node_to_attached_vec[n];
    }
    const OttIdSet & get_intruding_taxa() const {
        return intruding_taxa;
    }
    void set_intruding_taxa(const OttIdSet & ois) {
        intruding_taxa = ois;
    }
    private:
    long ott_id = 0L;
    const N * mrca = nullptr;
    // for each attachment point, we store which subtrees for this taxon attach there.
    attachment_to_node_vec attach_node_to_attached_vec;
    OttIdSet intruding_taxa;
};


template<typename N>
inline void writeLostTaxa(std::ostream & out,
                          const std::map<OttId, LostTaxonDetails<N> > & ltm, 
                          const std::map<OttId, OttIdSet> & synmap) {
    using json = nlohmann::json;
    json lostTaxa;
    for (auto ltmEl : ltm) {
        const auto & ottId = ltmEl.first;
        const auto & ltl = ltmEl.second; 
        json lostTaxonLocation;
        lostTaxonLocation["mrca"] = ltl.get_mrca()->get_name();
        const auto & attachmentPoints = ltl.get_attach_node_to_attached_vec();
        json apjson;
        for (auto attachPoint: attachmentPoints) {
            const auto & parent = attachPoint.first;
            const auto & childVec = attachPoint.second;
            json childVecJSON = json::array();
            for (auto c : childVec) {
                childVecJSON.push_back(c->get_name());
            }
            apjson[parent->get_name()] = childVecJSON;
        }
        lostTaxonLocation["attachment_points"] = apjson;
        const auto & intruding_taxa = ltl.get_intruding_taxa();
        if (!intruding_taxa.empty()) {
            json itj;
            for (auto i : intruding_taxa) {
                itj.push_back(i);
            }
            lostTaxonLocation["intruding_taxa"] = itj;
        }
        lostTaxa["ott" + std::to_string(ottId)] = lostTaxonLocation;
    }
    json synTaxa;
    for (auto synEl : synmap) {
        const auto & ottId = synEl.first;
        const auto & synset = synEl.second;
        json synj;
        for (auto jsEl : synset) {
            synj.push_back("ott" + std::to_string(jsEl));
        }
        synTaxa["ott" + std::to_string(ottId)] = synj;
    }
    json document;
    document["non_monophyletic_taxa"] = lostTaxa;
    document["taxa_matching_multiple_ott_ids"] = synTaxa;
    out << document.dump(1) << std::endl;
}


// Walk tipward from the @mrca until we find subtrees whose
// leaves are purely descendants of @taxon.  Record these
// subtrees as the attachment points in @ltl.
template <typename N>
inline void register_attachment_points(const N * taxon,
                                     const N * mrca,
                                     LostTaxonDetails<N> & ltl) {
    assert(taxon);
    assert(mrca);
    std::set<const N *> visited;
    std::stack<const N *> toDealWith;
    LOG(DEBUG) << "register_attachment_points for " << taxon->get_ott_id();
    const auto & taxDes = taxon->get_data().des_ids;
    assert(taxDes.size() > 1);
    assert(is_proper_subset(taxDes, mrca->get_data().des_ids));
    const N * currSolnNd = mrca;
    // the root of the solution, must have all of the constituent taxa
    while (true) {
        if (visited.count(currSolnNd) > 0) {
            if (toDealWith.empty()) {
                return ;
            }
            currSolnNd = toDealWith.top();
            toDealWith.pop();
            assert (visited.count(currSolnNd) == 0);
        } else {
            visited.insert(currSolnNd);
            const auto & solDes = currSolnNd->get_data().des_ids;
            assert(solDes != taxDes); // we should not call this if a des of mrca = taxon
            for (auto c : iter_child(*currSolnNd)) {
                const auto & cdesId = c->get_data().des_ids;
                if (!are_disjoint(cdesId, taxDes)) {
                    if (is_proper_subset(cdesId, taxDes)) {
                        ltl.addAttachmentPoint(currSolnNd, c);
                    } else {
                        toDealWith.push(c);
                    }
                }
            }
        }
    }
}

/// In our current impl. of the unpruner, the broken incertae sedis taxa
//    are hard to register in a lost taxon map during the unpruning.
//  `ott_id` is the ID of the "lost" taxon
//  `des_ids_for_taxon` should hold the IDs of all of its descendants.
//  `ott_to_node` is an ID to node map used to navigate the tree that breaks taxon
//  `lost_taxa_map` is the mapping to be modified. Overrides any entry in that 
//    map for `ott_id`
//  Note is OK for `des_ids_for_taxon` to include higher taxa IDs which are not
//    found in ott_to_node
template<typename N>
inline void fix_lost_taxon_map(OttId ott_id,
                               const OttIdSet & des_ids_for_taxon,
                               const std::map<OttId, N *> & ott_to_node,
                               std::map<OttId, LostTaxonDetails<N> > & lost_taxa_map) {
    auto tips = get_values_for_found_keys(ott_to_node, des_ids_for_taxon);
    std::set<N *> visited; // will hold pointers to all nodes in the induced tree
    const N * mrca = find_mrca_via_traversing(tips, &visited);
    // will hold pointers to all of the nodes in the summary tree that are
    //    in the induced tree and only have children that are also in the induced
    //    tree 
    std::set<const N *> no_interloper; 
    for (auto n: visited) {
        if (tips.count(n) == 0) {
            bool has_interloper = false;
            for (auto c : iter_child(*n)) {
                if (visited.count(c) == 0) {
                    has_interloper = true;
                    break;
                }
            }
            if (!has_interloper) {
                no_interloper.insert(n);
            }
        }
    }
    // Walking back from each tip, we compile a set of nodes that have children
    //    in and out of the taxon. These will be the MRCA or nodes not in no_interloper
    std::set<const N *> attachment_set;
    for (auto n : tips) {
        if (n == mrca) {
            attachment_set.insert(mrca);
            continue;
        }
        auto p = n->get_parent();
        while (p != mrca && no_interloper.count(p) > 0) {
            p = p->get_parent();
        }
        attachment_set.insert(p);
    }
    // compile a list of OTT IDs of intruding taxa, by getting the OTT Ids from
    // each attachment child that is not wholly assigned to this taxon
    // this traversal can include too many IDs we'll fix that below...
    OttIdSet intruding_taxa_plus;
    for (auto nd : attachment_set) {
        for (auto c : iter_child_const(*nd)) {
            if (no_interloper.count(c) == 0) {
                accumulate_closest_ott_id_for_subtree(c, intruding_taxa_plus);
            }
        }
    }
    const auto intruding_taxa = set_difference_as_set(intruding_taxa_plus, des_ids_for_taxon);
    // now we overwrite the existing entry with a clean LostTaxonDetails object, and fill it.
    lost_taxa_map[ott_id] = LostTaxonDetails<N>(ott_id, mrca);
    LostTaxonDetails<N> & ltl = lost_taxa_map[ott_id];
    ltl.set_intruding_taxa(intruding_taxa);
    // Fill in the attachments vector with all of the children that are tips of this taxon
    //    or "no_interloper" elements
    auto & attachment = ltl.get_attach_node_to_attached_vec();
    attachment.clear();
    for (auto n: attachment_set) {
        auto & att_vec = attachment[n];
        for (auto c : iter_child(*n)) {
            if (tips.count(const_cast<N *>(c)) > 0 || no_interloper.count(c) > 0) {
                att_vec.insert(c);
            }
        }
    }
}


} //namespace otc
#endif
