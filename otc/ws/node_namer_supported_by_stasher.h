#ifndef OTC_NODE_NAMER_STASHER_H
#define OTC_NODE_NAMER_STASHER_H
#include "otc/ws/tolws.h"
#include "otc/ws/trees_to_serve.h"


namespace otc {



inline std::string node_id_for_summary_tree_node(const SumTreeNode_t & nd) {
    //return nd.get_name();
    // code below from when we were trying calling clear_names_for_nodes_with_ids to clear
    //    the name strings, but I that saved trivial memory
    return (nd.has_ott_id() ? ott_id_to_idstr(nd.get_ott_id()) : nd.get_name());
}

class NodeNamerSupportedByStasher {
    public:
    mutable std::set<const std::string *> study_id_set;
    NodeNameStyle nns;
    const RichTaxonomy & taxonomy;
    const TreesToServe & tts;
    NodeNamerSupportedByStasher(NodeNameStyle in_nns, const RichTaxonomy &tax, const TreesToServe & tts_arg)
        :nns(in_nns),
        taxonomy(tax),
        tts(tts_arg) {
    }

    std::string operator()(const SumTreeNode_t *nd) const {
        const SumTreeNodeData & d = nd->get_data();
#           if defined(JOINT_MAPPING_VEC)
            for (auto & el : d.source_edge_mappings) {
                if (el.first == SourceEdgeMappingType::SUPPORTED_BY_MAPPING) {
                    study_id_set.insert(tts.decode_study_node_id_index(el.second).first);
                }
            }
#           else
            for (auto & sni : d.supported_by) {
                study_id_set.insert(tts.decode_study_node_id_index(sni).first);
            }
#           endif
        if (nns != NodeNameStyle::NNS_ID_ONLY && nd->has_ott_id()) {
            const auto * tr = taxonomy.included_taxon_from_id(nd->get_ott_id());
            if (tr == nullptr) {
                throw OTCError() << "OTT Id " << nd->get_ott_id() << " in namer not found in taxonomy! Please report this bug";
            }
            std::string taxon_name = get_taxon_unique_name(*tr);
            if (nns == NodeNameStyle::NNS_NAME_AND_ID) {
                std::string ret;
                const std::string id_str = node_id_for_summary_tree_node(*nd);
                ret.reserve(taxon_name.length() + 1 + id_str.length());
                ret = taxon_name;
                ret += ' ';
                ret += id_str;
                return ret;    
            } else {
                return taxon_name;
            }
        }
        return node_id_for_summary_tree_node(*nd);
    }

    std::string operator()(const RTRichTaxNode *nd) const {
        assert(nd != nullptr);
        if (nns == NodeNameStyle::NNS_ID_ONLY) {
            return ott_id_to_idstr(nd->get_ott_id());
        }
        const std::string & taxon_name = get_taxon_unique_name(*nd);
        if (nns == NodeNameStyle::NNS_NAME_AND_ID) {
            std::string ret;
            std::string id_str = ott_id_to_idstr(nd->get_ott_id());
            ret.reserve(taxon_name.length() + 1 + id_str.length());
            ret = taxon_name;
            ret += ' ';
            ret += id_str;
            return ret;    
        }
        return taxon_name;
    }
};


} // namespace otc
#endif
