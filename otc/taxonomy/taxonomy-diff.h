#ifndef OTC_TAXONOMY_DIFF_H
#define OTC_TAXONOMY_DIFF_H

#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <bitset>
#include <regex>


#include "json.hpp"

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/diff_maker.h"
#include "otc/config_file.h"

namespace otc {
template<typename T>
bool fill_name_id_maps(const T & tree_data,
                      std::unordered_map<OttId, std::string_view>& id2name) {
    for (auto id_nd_pair : tree_data.id_to_node) {
        const auto & nd_ptr = id_nd_pair.second;
        const auto tax_id = id_nd_pair.first;
        if (nd_ptr == nullptr) {
            throw OTCError() << "Unexpected nullptr in id_to_node";
        }
        id2name[tax_id] = nd_ptr->get_name();
    }
    return true;
}

inline
OttIdSet find_ids_with_same_names(const std::unordered_map<OttId, std::string_view> & old_id2name, 
                                  const std::unordered_map<OttId, std::string_view> & new_id2name) {
    OttIdSet same_id_name;
    for (auto id_name_pair : old_id2name) {
        OttId tax_id = id_name_pair.first;
        auto new_it = new_id2name.find(tax_id);
        if (new_it != new_id2name.end()) {
            if (new_it->second == id_name_pair.second) {
                same_id_name.insert(tax_id);
            }
        }
    }
    return same_id_name;
}

inline
std::unordered_map<const RTRichTaxNode *, OttIdSet> fill_term_des_id_set(const RichTaxTree & tree,
                                const OttIdSet & relevant_ids,
                                OttIdSet & terms_added) {
    std::unordered_map<const RTRichTaxNode *, OttIdSet> nd2idset;
    for (auto nd : iter_post_const(tree)) {
        const OttId tax_id = nd->get_ott_id();
        const bool is_relevant = contains(relevant_ids, tax_id);
        OttIdSet & dest = nd2idset[nd];
        // add a tip if it is not the ancestor of another relevant ID
        if (!nd->is_tip()) {
            for (auto child : iter_child_const(*nd)) {
                const auto & child_set = nd2idset.at(child);
                dest.insert(child_set.begin(), child_set.end());
            }
        }
        if (is_relevant && dest.empty()) {
            dest.insert(tax_id);
            terms_added.insert(tax_id);
        }
    }
    return nd2idset;
}

inline
void partitionTaxonByTypeOfType(const RichTaxTree & tree,
                                OttIdSet & clade_ids,
                                OttIdSet & specimen_based) {
    std::vector<const RTRichTaxNode *> unk;
    for (auto nd : iter_pre_const(tree)) {
        auto & nd_data = nd->get_data();
        auto & rank = nd_data.rank;
        bool is_clade = false;
        const auto tax_id = nd->get_ott_id();
        assert(!contains(clade_ids, tax_id));
        assert(!contains(specimen_based, tax_id));
        if (rank == TaxonomicRank::RANK_NO_RANK && !nd->is_tip()) {
            unk.push_back(nd);
            continue;
        }
        if (rank < TaxonomicRank::RANK_SPECIES) {
            is_clade = true;
        }
        const auto & par = nd->get_parent();
        if (par) {
            if (is_clade) {
                if (contains(specimen_based, par->get_ott_id())) {
                    throw OTCError() << "Taxon with ID " << tax_id << " is higher taxon, but parent (" << par->get_ott_id() << ") is not.";
                }
            }
        }
        if (is_clade) {
            clade_ids.insert(tax_id);
        } else {
            specimen_based.insert(tax_id);
        }
    }
    while (unk.size() > 0) {
        std::vector<const RTRichTaxNode *> current;
        current.clear();
        std::swap(unk, current);
        const auto before_size = current.size();
        for (auto nd : current) {
            const auto tax_id = nd->get_ott_id();
            bool is_below_sp = false;
            for (auto anc : iter_anc_const(*nd)) {
                auto anc_id = anc->get_ott_id();
                if (contains(specimen_based, anc_id)) {
                    // LOG(DEBUG) << "Anc " << anc_id << " causing " << tax_id << " to be specimen_based.";
                    is_below_sp = true;
                    break;
                }
                if (contains(clade_ids, anc_id)) {
                    break;
                }
            }
            if (is_below_sp) {
                specimen_based.insert(tax_id);
                continue;
            }
            bool is_above_higher = false;
            for (auto c : iter_child_const(*nd)) {
                auto c_id = c->get_ott_id();
                if (contains(clade_ids, c_id)) {
                    // LOG(DEBUG) << "child " << c_id << " causing " << tax_id << " to be clade.";
                    is_above_higher = true;
                    break;
                }
            }
            if (is_above_higher) {
                clade_ids.insert(tax_id);
            } else {
                unk.push_back(nd);
            }
        }
        const auto after_size = unk.size();
        if (before_size == after_size) {
            // creeping from edges has stopped settling taxa
            //  call all of the rest as higher taxa
            for (auto nd : current) {
                clade_ids.insert(nd->get_ott_id());
            }
            unk.clear();
        }
        
    }
}

            
inline
const RTRichTaxNode * find_specimen_based_root(const RTRichTaxNode *nd,
                                               const OttIdSet & spec_based_ids,
                                               OttIdSet & seen) {
    const RTRichTaxNode * par = nd->get_parent();
    const RTRichTaxNode * ret = nullptr;
    while (true) {
        if (!par) {
            ret = nd;
            break;
        }
        if (contains(spec_based_ids, par->get_ott_id())) {
            nd = par;
        } else {
            ret = nd;
            break;
        }
        par = nd->get_parent();
    }
    if (ret->is_tip()) {
        seen.insert(ret->get_ott_id());
    } else {
        for (auto nnd : iter_pre_n_const(ret)) {
            auto tax_id = nnd->get_ott_id();
            if (!contains(spec_based_ids, tax_id)) {
                throw OTCError() << "taxon " << nnd->get_name() << " (" << tax_id << ") not specimen_based, but in spec_based clade";
            }
            seen.insert(tax_id);
        }
    }
    return ret;
}

inline
std::unordered_map<OttId, const RTRichTaxNode *> find_all_specimen_based_roots(const std::unordered_map<OttId, const RTRichTaxNode *> & id2nd,
                                      const OttIdSet & specimen_based_ids,
                                      OttIdSet & seen) {
    std::unordered_map<OttId, const RTRichTaxNode *> sp_root;
    for (auto ott_id : specimen_based_ids) {
        if (contains(seen, ott_id)) {
            continue;
        }
        auto nd = id2nd.at(ott_id);
        auto root_ptr = find_specimen_based_root(nd, specimen_based_ids, seen);
        assert(root_ptr != nullptr);
        sp_root[root_ptr->get_ott_id()] = root_ptr;
    }
    return sp_root;
}


enum AlphaEditOp {
    NO_CHANGE = 0,
    CHANGED_ID = 1,
    CHANGED_NAME = 2,
    DELETED_SYN = 3,
    ADDED_SYN = 4,
    CHANGED_RANK = 5,
    DELETE_TAXON = 6,
    ADD_TAXON = 7,
    CHANGED_FLAGS = 8
};

const std::vector<std::string> aeo2str = {"no change", "change id", "change name",
                                "delete synonym", "add synonym", 
                                "change rank", "delete taxon", "add taxon",
                                "change flags"};

const std::map<std::string, AlphaEditOp> str2aeo = {
    {"no change", AlphaEditOp::NO_CHANGE},
    {"change id", AlphaEditOp::CHANGED_ID},
    {"change name", AlphaEditOp::CHANGED_NAME},
    {"delete synonym", AlphaEditOp::DELETED_SYN},
    {"add synonym", AlphaEditOp::ADDED_SYN},
    {"change rank", AlphaEditOp::CHANGED_RANK},
    {"delete taxon", AlphaEditOp::DELETE_TAXON},
    {"add taxon", AlphaEditOp::ADD_TAXON},
    {"change flags", AlphaEditOp::CHANGED_FLAGS}
};

class AlphaEdit {
    public:
        AlphaEdit()
            :operation(AlphaEditOp::NO_CHANGE),
            first_rank(TaxonomicRank::RANK_NO_RANK),
            second_rank(TaxonomicRank::RANK_NO_RANK) {
        }

        AlphaEditOp operation;
        std::string first_str, second_str;
        OttId first_id, second_id;
        TaxonomicRank first_rank, second_rank;
        tax_flags first_flags, second_flags;

        void add_to_json_array(nlohmann::json & jarr) const {
            if (operation == AlphaEditOp::NO_CHANGE) {
                return;
            } 
            nlohmann::json el;
            el["operation"] = aeo2str[operation];
            if (operation == AlphaEditOp::CHANGED_ID) {
                el["from"] = first_id;
                el["to"] = second_id;
            } else {
                el["taxon_id"] = first_id;
                if (operation == AlphaEditOp::CHANGED_NAME) {
                    el["from"] = first_str;
                    el["to"] = second_str;
                } else if (operation == AlphaEditOp::DELETED_SYN || operation == AlphaEditOp::ADDED_SYN) {
                    el["synonym"] = first_str;
                    if (!second_str.empty()) {
                        el["type"] = second_str;
                    }
                } else if (operation == AlphaEditOp::DELETE_TAXON) {
                    if (!first_str.empty()) {
                        el["name"] = first_str;
                    }
                } else if (operation == AlphaEditOp::CHANGED_RANK) {
                    el["from"] = rank_enum_to_name.at(first_rank);
                    el["to"] = rank_enum_to_name.at(second_rank);
                } else if (operation == AlphaEditOp::CHANGED_FLAGS) {
                    el["from"] = flags_to_string(first_flags);
                    el["to"] = flags_to_string(second_flags);
                } else {
                    assert(operation == AlphaEditOp::ADD_TAXON);
                    el["name"] = first_str;
                    el["rank"] = rank_enum_to_name.at(first_rank);
                    auto fstr = flags_to_string(first_flags);
                    if (!fstr.empty()) {
                        el["flags"] = fstr;
                    }
                }
            }
            jarr.push_back(el);
        }
};

enum AlphaGroupEditOp {
    NO_GR_CHANGE = 0,
    GR_CHANGED_ID = 1,
    GR_CHANGED_NAME = 2,
    ADD_TAXA = 3,
    // DEL_TAXA = 4,
    // ADD_DEL_TAXA = 5,
    DELETED_GROUPING = 6,
    NEW_GROUPING = 7,
    GR_CHANGED_RANK = 8,
    GR_CHANGED_FLAGS = 9
};

const std::vector<std::string> ageo2str = {"no change", "change id", "change name",
                                "add taxa", "delete taxa", "add+delete taxa",
                                "deleted grouping", "new grouping",
                                "change rank", "change flags"};

const std::map<std::string, AlphaGroupEditOp> str2ageo = {
    {"no change", AlphaGroupEditOp::NO_GR_CHANGE},
    {"change id", AlphaGroupEditOp::GR_CHANGED_ID},
    {"change name", AlphaGroupEditOp::GR_CHANGED_NAME},
    {"add taxa", AlphaGroupEditOp::ADD_TAXA},
    // {"delete taxa", AlphaGroupEditOp::DEL_TAXA},
    // {"add+delete taxa", AlphaGroupEditOp::ADD_DEL_TAXA},
    {"deleted grouping", AlphaGroupEditOp::DELETED_GROUPING},
    {"new grouping", AlphaGroupEditOp::NEW_GROUPING},
    {"change rank", AlphaGroupEditOp::GR_CHANGED_RANK},
    {"change flags", AlphaGroupEditOp::GR_CHANGED_FLAGS}

};

class AlphaGroupEdit {
    public:
        AlphaGroupEdit()
            :operation(AlphaGroupEditOp::NO_GR_CHANGE) {
        }

        AlphaGroupEditOp operation;
        std::string first_str, second_str;
        OttId first_id, second_id;
        OttIdSet newChildIds;
        //OttIdSet addedIds;
        //OttIdSet delIds;
        TaxonomicRank first_rank, second_rank;
        tax_flags first_flags, second_flags;

        void add_to_json_array(nlohmann::json & jarr) const {
            if (operation == AlphaGroupEditOp::NO_GR_CHANGE) {
                return;
            } 
            nlohmann::json el;
            el["operation"] = ageo2str[operation];
            if (operation == AlphaGroupEditOp::GR_CHANGED_ID) {
                el["from"] = first_id;
                el["to"] = second_id;
            } else {
                el["taxon_id"] = first_id;
                if (operation == AlphaGroupEditOp::GR_CHANGED_NAME) {
                    el["from"] = first_str;
                    el["to"] = second_str;
                } else if (operation == AlphaGroupEditOp::GR_CHANGED_RANK) {
                    el["from"] = rank_enum_to_name.at(first_rank);
                    el["to"] = rank_enum_to_name.at(second_rank);
                } else if (operation == AlphaGroupEditOp::GR_CHANGED_FLAGS) {
                    el["from"] = flags_to_string(first_flags);
                    el["to"] = flags_to_string(second_flags);
                } else if (operation == AlphaGroupEditOp::ADD_TAXA) {
                    nlohmann::json added = nlohmann::json::array();
                    for (auto oid : newChildIds) {
                        added.push_back(oid);
                    }
                    el["added"] = added;
                // } else if (operation == AlphaGroupEditOp::DEL_TAXA) {
                //     nlohmann::json deleted = nlohmann::json::array();
                //     for (auto oid : delIds) {
                //         deleted.push_back(oid);
                //     }
                //     el["deleted"] = deleted;
                // } else if (operation == AlphaGroupEditOp::ADD_DEL_TAXA) {
                //     nlohmann::json added = nlohmann::json::array();
                //     for (auto oid : addedIds) {
                //         added.push_back(oid);
                //     }
                //     el["added"] = added;
                //     nlohmann::json deleted = nlohmann::json::array();
                //     for (auto oid : delIds) {
                //         deleted.push_back(oid);
                //     }
                //     el["deleted"] = deleted;
                } else if (operation == AlphaGroupEditOp::NEW_GROUPING) {
                    nlohmann::json added = nlohmann::json::array();
                    for (auto oid : newChildIds) {
                        added.push_back(oid);
                    }
                    el["added"] = added;
                    el["name"] = first_str;
                    el["rank"] = rank_enum_to_name.at(first_rank);
                    auto fstr = flags_to_string(first_flags);
                    if (!fstr.empty()) {
                        el["flags"] = fstr;
                    }
                } else if (operation == AlphaGroupEditOp::DELETED_GROUPING) {
                    // no-op
                }
            }
            jarr.push_back(el);
        }
};

class Grouping {
    public:
    Grouping() : paired(false){}

    std::string name;
    OttId tax_id;
    OttIdSet shared_ids; // always put in terms of "new" ids
    OttIdSet new_ids;    // new taxa in new group;
    bool paired;
    const RTRichTaxNode * node;
};

} // namespace


#endif