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
///////////////////////////////////////////////////////////////////////////////
// LostTaxonDetails holds the information about one of the broken taxa.
// 
// It stores info for a taxon X which not in the solution.
//   1. get_mrca() returns a pointer to the MRCA of the constituents of the taxon
//   2. get_attach_node_to_attached_set() returns a mapping between nodes
//        that are in the subtree of the MRCA and which contain a mixture
//        of some children (not listed) that are not purely of this taxon and
//        some children (in the vector that is the value) that contain only 
//        descendants that are labelled as being a part of taxon X.
//   3. get_intruding_taxa() returns a set of OTT Ids that are excluded from X
//        by the defition of X, but are descendants of the MRCA. Note that if
//        a higher level taxon is an "intruder" then only that taxon's ID (not
//        all of its descendants' IDs) will be listed.
//   4. get_allowed_intruding_taxa() lists intruders that are not excluded because
//        they are incertae sedis taxa (or descendants of incertae sedis). This
//        list is not exhaustive. If an incertae sedis taxon is placed inside a higher
//        taxon that is a descendant of taxon X, then the i.s. taxon will not be
//        included. So this only refers to the i.s. taxa that attach within the scope
//        of the attachment points of this taxon.
// 
// An attachment point is:
//      1. an ancestor of at least one member of the taxon X,
//      2. an ancestor of at least one taxon that is not contained in the taxon X.
//      3. has at least one child (an attached node) that consists entirely
//          of member of taxon X
// The last condition means that in:
//     ((X1,O1)N1,(X2,O2)N2)N3
// N1 and N2 are attachement points for taxon X, but N3 is not.
template<typename N>
class LostTaxonDetails {
    using node_set = std::set<const N *>;
    using attachment_to_node_set = std::map<const N *, node_set>;
    public:
    LostTaxonDetails& operator=(LostTaxonDetails&& tr) = default;
    LostTaxonDetails(LostTaxonDetails&& tr) = default;
    
    LostTaxonDetails(OttId oid,
                     const OttIdSet & des_ids_for_taxon,
                     const OttIdSet & non_excluded_ids_for_taxon,
                     const std::map<OttId, N *> & ott_to_node);
    const N * get_mrca() const {
        return mrca;
    }
    const attachment_to_node_set & get_attach_node_to_attached_set() const {
        return attach_node_to_attached_set;
    }
    const node_set & get_attached_set(const N * n) const {
        return attach_node_to_attached_set.at(n);
    }
    const OttIdSet & get_intruding_taxa() const {
        return intruding_taxa;
    }
    const OttIdSet & get_allowed_intruding_taxa() const {
        return allowed_intruding_taxa;
    }
    node_set & get_attached_set(const N * n) {
        return attach_node_to_attached_set[n];
    }
    private:
    void add_child_to_attachment_point(const N * attach, const N * child) {
        attach_node_to_attached_set[attach].insert(child);
    }
    attachment_to_node_set & get_attach_node_to_attached_set() {
        return attach_node_to_attached_set;
    }
    long ott_id = 0L;
    const N * mrca = nullptr;
    // for each attachment point, we store which subtrees for this taxon attach there.
    attachment_to_node_set attach_node_to_attached_set;
    OttIdSet intruding_taxa;
    OttIdSet allowed_intruding_taxa; // incertae sedis taxa that can invade without breaking.
};


template<typename N>
inline void write_lost_taxa(std::ostream & out,
                            const std::map<OttId, LostTaxonDetails<N> > & ltm, 
                            const std::map<OttId, OttIdSet> & synmap) {
    using json = nlohmann::json;
    json lostTaxa;
    for (auto ltm_it = ltm.begin(); ltm_it != ltm.end(); ++ltm_it) {
        const auto & ottId = ltm_it->first;
        const auto & ltl = ltm_it->second; 
        json lostTaxonLocation;
        lostTaxonLocation["mrca"] = ltl.get_mrca()->get_name();
        const auto & attachmentPoints = ltl.get_attach_node_to_attached_set();
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
                itj.push_back("ott" + std::to_string(i));
            }
            lostTaxonLocation["intruding_taxa"] = itj;
        }
        const auto & allowed_intruding_taxa = ltl.get_allowed_intruding_taxa();
        if (!allowed_intruding_taxa.empty()) {
            json itj;
            for (auto i : allowed_intruding_taxa) {
                itj.push_back("ott" + std::to_string(i));
            }
            lostTaxonLocation["nonexcluded_intruding_taxa"] = itj;
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

/// In our current impl. of the unpruner, the broken incertae sedis taxa
//    are hard to register in a lost taxon map during the unpruning.
// This function fills in a LostTaxonDetails struct for a taxon without requiring a
//    des_ids field for the tree.
//  `ott_id` is the ID of the "lost" taxon
//  `des_ids_for_taxon` should hold the IDs of all of its descendants.
//  `ott_to_node` is an ID to node map used to navigate the tree that breaks taxon
//  `lost_taxa_map` is the mapping to be modified. Overrides any entry in that 
//    map for `ott_id`
//  Note is OK for `des_ids_for_taxon` to include higher taxa IDs which are not
//    found in ott_to_node
template<typename N>
inline LostTaxonDetails<N>::LostTaxonDetails(OttId oid,
                                             const OttIdSet & des_ids_for_taxon,
                                             const OttIdSet & non_excluded_ids_for_taxon,
                                             const std::map<OttId, N *> & ott_to_node)
     :ott_id(oid),
      mrca(nullptr) {
    auto tips = get_values_for_found_keys(ott_to_node, des_ids_for_taxon);
    std::set<N *> visited; // will hold pointers to all nodes in the induced tree
    this->mrca = find_mrca_via_traversing(tips, &visited);
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
    auto intruding_taxa_uncl = set_difference_as_set(intruding_taxa_plus, des_ids_for_taxon);
    this->intruding_taxa = set_difference_as_set(intruding_taxa_uncl, non_excluded_ids_for_taxon);
    this->allowed_intruding_taxa = set_intersection_as_set(intruding_taxa_uncl, non_excluded_ids_for_taxon);
    // Fill in the attachments vector with all of the children that are tips of this taxon
    //    or in the no_interloper set
    auto & attachment = this->get_attach_node_to_attached_set();
    attachment.clear();
    for (auto n: attachment_set) {
        auto & att_set = attachment[n];
        for (auto c : iter_child(*n)) {
            if (tips.count(const_cast<N *>(c)) > 0 || no_interloper.count(c) > 0) {
                att_set.insert(c);
            }
        }
    }
}


} //namespace otc
#endif
