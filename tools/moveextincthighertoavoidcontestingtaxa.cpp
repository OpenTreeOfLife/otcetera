#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include "json.hpp"
#include <sstream>
using namespace otc;
using json = nlohmann::json;

bool read_json_array(json & j, std::set<OttId> & set_ints);


bool read_json_array(json & j, std::set<OttId> & set_ints) {
    if (j.type() != json::value_t::array) {
        throw OTCError("Expecting array in json");
    }
    unsigned index = 0;
    for (auto i : j) {
        auto itype = i.type();
        if (itype != json::value_t::number_integer && itype != json::value_t::number_unsigned) {
            LOG(ERROR) << "Could not parse all elements of input list of as integer:\n" << i ;
        return false;
        }
        unsigned val = i.get<OttId>();
        //std::cerr << "Element[" << index++ << "] = " << val << '\n';
        set_ints.insert(val);
    }
    return true;
}


struct MoveExtinctHigherState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    using NeedsMoveTipmostTaxonPair = std::pair<bool, const NodeWithSplits *>;
    using TreeIDToPlacement = std::map<std::string, NeedsMoveTipmostTaxonPair>;
    using ExtinctTaxonToPlacementSummary = std::map<OttId, TreeIDToPlacement>;
    using ConstNdPtr = const TreeMappedWithSplits::node_type *;
    using Phylo2Taxo = std::map<ConstNdPtr, ConstNdPtr>;
    
    using InducedParAndIds = std::pair<ConstNdPtr, OttIdSet>;
    // used to map taxonomy node (x) to induced tree parent node (y) and y's des_ids restricted to the set of tips in the phylo
    using TaxoNd2IndParIDset = std::map<ConstNdPtr,  InducedParAndIds> ;
    
    using TaxoNd2RestrictedIDSet = std::map<ConstNdPtr, OttIdSet>;
    using ConstNdPtrSet = std::set<ConstNdPtr>;
    using ConstNdPtrPair = std::pair<ConstNdPtr, ConstNdPtr> ;
    using DeepestAttachmentPointAndTreeIDs = std::pair<const NodeWithSplits *, std::set<std::string> >;
    using ExtinctToDeepest = std::map<const NodeWithSplits *, DeepestAttachmentPointAndTreeIDs>;

    int numErrors;
    bool useStdOut;
    bool useCmdLine;
    std::string extinctInputJSONFilepath;
    std::string extinctInputJSONInline;
    std::string outputJSONLogFile;
    std::istream * extinctInStream;
    std::ostream * jsonOutStream;
    std::ofstream jsonOutFileIfUsed;
    std::set<OttId> extinctIDSet;
    std::set<OttId> incertSedisIDSet;
    std::size_t phylogeniesProcessed = 0;
    std::vector<std::string> namesOfTreesWithoutExtinct;
    ExtinctTaxonToPlacementSummary unjoinedExtinctTaxonToPlacementSummary; 
    ExtinctToDeepest joinedExtinctTaxonToPlacementSummary;

    const NodeWithSplits * find_shallowest_fossil_only(const NodeWithSplits *q) {
        const NodeWithSplits * curr = q;
        if (!contains(extinctIDSet, q->get_ott_id())) {
            throw OTCError() << "Called find_shallowest_fossil_only with non extinct node, with OTT ID = " << q->get_ott_id();
        }
        for (;;) {
            const NodeWithSplits * next_nd = curr->get_parent();
            if (next_nd == nullptr || !contains(extinctIDSet, next_nd->get_ott_id())) {
                return curr;
            }
            curr = next_nd;
        }
    }

    DeepestAttachmentPointAndTreeIDs find_deepest_attachment(const TreeIDToPlacement & edits) {
        std::set<std::string> tree_names;
        const NodeWithSplits * deepest_attach_nd = nullptr;
        int lowest_depth = 0;
        for (auto tree_pref_pair: edits) {
            auto needs_move_taxon = tree_pref_pair.second;
            if (needs_move_taxon.first) {
                bool lowest = false;
                auto attach_nd = needs_move_taxon.second;
                const RTSplits & an_data = attach_nd->get_data();
                auto an_depth = an_data.depth;
                if (tree_names.empty() || an_depth < lowest_depth) {
                    lowest = true;
                    lowest_depth = an_depth;
                    tree_names.clear();
                    deepest_attach_nd = attach_nd;
                }
                if (lowest_depth == an_depth) {
                    tree_names.insert(tree_pref_pair.first);
                }
            }
        }
        return DeepestAttachmentPointAndTreeIDs{deepest_attach_nd, tree_names};
    }

    void joinPlacements() {
        const RTreeOttIDMapping<RTSplits> & taxdata = taxonomy->get_data();
        // unjoinedButDeepest will map all nodes that need to be moved to their new 
        //    attachment points (and keep list of tree names)
        ExtinctToDeepest unjoinedButDeepest;
        for (auto tax_id_to_edit : unjoinedExtinctTaxonToPlacementSummary) {
            const auto edits = tax_id_to_edit.second;
            auto deepest_and_prov = find_deepest_attachment(edits);
            if (deepest_and_prov.first == nullptr) {
                continue;
            }
            const auto moving_ott_id = tax_id_to_edit.first;
            const auto taxon_ptr = taxdata.ott_id_to_node.at(moving_ott_id);
            unjoinedButDeepest[taxon_ptr] = deepest_and_prov;
        }

        // keys of obs_to_fossil_anc will be the most inclusive fossil ancestors that need to move
        // values will be the set of keys in unjoinedButDeepest
        std::map<const NodeWithSplits *, std::set<const NodeWithSplits *> > obs_to_fossil_anc;
        for (auto fnd : unjoinedButDeepest) {
            auto obs_node = fnd.first;
            auto af_node = find_shallowest_fossil_only(obs_node);
            obs_to_fossil_anc[af_node].insert(obs_node);
        }

        // Now for each fossil taxon that is moving, choose the deepest attachament point
        joinedExtinctTaxonToPlacementSummary.clear();
        for (auto otfa_it : obs_to_fossil_anc) {
            auto moving_nd = otfa_it.first;
            const std::set<const NodeWithSplits *> & obs_set = otfa_it.second;
            int ld = 0;
            DeepestAttachmentPointAndTreeIDs dapati;
            dapati.first = nullptr;
            for (auto obs : obs_set) {
                const auto & deepest_for_obs = unjoinedButDeepest.at(obs);
                auto attach_nd = deepest_for_obs.first;
                const RTSplits & an_data = attach_nd->get_data();
                auto an_depth = an_data.depth;
                if (dapati.first == nullptr || an_depth < ld) {
                    dapati = deepest_for_obs;
                    ld = an_depth;
                } else if (an_depth < ld) {
                    const std::set<std::string> & new_tree_names = deepest_for_obs.second;
                    dapati.second.insert(begin(new_tree_names), end(new_tree_names));
                }
            }
            joinedExtinctTaxonToPlacementSummary[moving_nd] = dapati;
        }
        unjoinedExtinctTaxonToPlacementSummary.clear();
    }
  
    /*

            std::set<std::string> 
        json empty;
        if (!(tree_names.empty()) {
            json jtree_names = json::array();
            for (auto i : tree_names) {
                jtree_names.push_back(i);
            }
            return AttachNameJSONPair(attach_nd_id, jtree_names);
        }
        return AttachNameJSONPair(std::string(), empty);
    }*/

    virtual ~MoveExtinctHigherState(){}
    MoveExtinctHigherState()
        :numErrors(0),
        useStdOut(true),
        useCmdLine(false) {
    }
    bool readExtinctIDsFromStream(OTCLI & otCLI, std::istream &inp) {
        json j;
        try {
            inp >> j;
        } catch (json::parse_error & x) {
            LOG(ERROR) << "Could not parse input JSON list of IDs.:\n" << x.what() << " at byte " << x.byte ;
            return false;
        }
        if (j.type() == json::value_t::object) {
            try {
                json earr = j.at("extinct");
                if (!read_json_array(earr, extinctIDSet)) {
                    return false;
                }
            } catch (json::out_of_range & ) {
            }
            try {
                auto isarr = j.at("incertae_sedis");
                if (!read_json_array(isarr, incertSedisIDSet)) {
                    return false;
                }
            } catch (json::out_of_range & ) {
            }
        }
        if (extinctIDSet.empty()) {
            LOG(ERROR) << "Expecting an object with an array of \"extinct\" OTT Ids as the JSON content.";
            return false;
        }
        
        return true;
    }

    bool readExtinctIDs(OTCLI & otCLI) {
        std::ifstream extinctJSONStreamLocal;
        std::istringstream extinctJSONStringWrapper;
        if (useCmdLine) {
            if (!extinctInputJSONFilepath.empty()) {
                LOG(ERROR) << "Cannot use command line JSON and input JSON file options at the same time!";
                return false;
            }
            extinctJSONStringWrapper.str(extinctInputJSONInline);
            extinctInStream = &extinctJSONStringWrapper;
        } else {
            if (extinctInputJSONFilepath.empty()) {
                LOG(ERROR) << "Must use either the command line JSON or input JSON file options!";
                return false;
            }
            extinctJSONStreamLocal.open(extinctInputJSONFilepath);
            if (!extinctJSONStreamLocal.good()) {
                LOG(ERROR) << "Could not open input extinct JSON file \"" << extinctInputJSONFilepath << "\"";
                return false;
            }
            extinctInStream = &extinctJSONStreamLocal;
        }
        if (useStdOut) {
            jsonOutStream = &std::cout;
        } else {
            jsonOutFileIfUsed.open(outputJSONLogFile);
            if (!jsonOutFileIfUsed.good()) {
                jsonOutFileIfUsed.close();
                LOG(ERROR) << "Could not open JSON output path at \"" << outputJSONLogFile << "\"";
                return false;
            }
            jsonOutStream = &jsonOutFileIfUsed;
        }
        return readExtinctIDsFromStream(otCLI, *extinctInStream);
    }

    bool pretree_read_hook(OTCLI & otCLI) override {
        return readExtinctIDs(otCLI);
    }
    
    bool process_taxonomy_tree(OTCLI & otCLI) override {
        bool r = TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::process_taxonomy_tree(otCLI);
        // we can ignore the internal node labels for the non-taxonomic trees
        otCLI.get_parsing_rules().set_ott_idForInternals = false;
        for (auto nd : iter_pre(*taxonomy)) {
            RTSplits & d = nd->get_data();
            auto p = nd->get_parent();
            if (p == nullptr) {
                d.depth = 0;
            } else {
                d.depth = 1 + p->get_data().depth;
            }
            //db_write_ott_id_set(nd->get_name().c_str(), d.des_ids);
            //std::cerr << "nd id = " << nd->get_ott_id() << " depth = " << d.depth << '\n';
        }
        return r;
    }

  
    bool summarize(OTCLI &otCLI) override {
        std::cerr << "unjoinedExtinctTaxonToPlacementSummary.size() = " << unjoinedExtinctTaxonToPlacementSummary.size() << '\n';
        joinPlacements();
        std::cerr << "unjoinedExtinctTaxonToPlacementSummary.size() after = " << unjoinedExtinctTaxonToPlacementSummary.size() << '\n';
        json document;
        json taxon_edits;
        std::cerr << "joinedExtinctTaxonToPlacementSummary.size() = " << joinedExtinctTaxonToPlacementSummary.size() << '\n';
        for (auto tax_id_to_edit : joinedExtinctTaxonToPlacementSummary) {
            const auto taxon_ptr = tax_id_to_edit.first;
            const DeepestAttachmentPointAndTreeIDs & edits = tax_id_to_edit.second;
            const auto & att_nd_id = edits.first;
            const auto & tree_set = edits.second;
            json new_attach;
            std::string par_key = "ott" + std::to_string(att_nd_id->get_ott_id());
            new_attach["new_parent"] = par_key;
            new_attach["according_to"] = tree_set;
            std::string key = "ott" + std::to_string(taxon_ptr->get_ott_id());
            taxon_edits[key] =  new_attach;
        }
        document["edits"] = taxon_edits;
        json jTreesWoExt = json::array(); 
        for (auto i : namesOfTreesWithoutExtinct) {
            jTreesWoExt.push_back(i);
        }
        document["trees_with_no_extinct_tips"] = jTreesWoExt;
        assert(jsonOutStream != nullptr);
        *jsonOutStream << document << std::endl;
        if (jsonOutStream == & jsonOutFileIfUsed) {
            jsonOutFileIfUsed.close();
            jsonOutStream = nullptr;  
        }
        return true;
    }
 

    bool process_source_tree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> treeup) override {
        std::size_t treeIndex = phylogeniesProcessed++;
        TreeMappedWithSplits * raw = treeup.get();
        raw->set_name(otCLI.currentFilename);
        ConstNdPtrSet extinct_tips;
        OttIdSet relExtinctIDs;
        OttIdSet relIncSedIDs;

        for (auto nd : iter_leaf(*raw)) {
            auto ottId = nd->get_ott_id();
            if (contains(extinctIDSet, ottId)) {
                relExtinctIDs.insert(ottId);
                extinct_tips.insert(nd);
                const RTSplits & d = nd->get_data();
                LOG(INFO) << "Tree " << treeIndex << " name=" <<  raw->get_name() << " has extinct tip with ID=" << ottId << ".\n";
            } else if (contains(incertSedisIDSet, ottId)) {
                relIncSedIDs.insert(ottId);
            }
        }
        if (extinct_tips.empty()) {
            LOG(INFO) << "Tree " << treeIndex << " name=" <<  raw->get_name() << " contains no extinct tips. flushing...\n";
            namesOfTreesWithoutExtinct.push_back(raw->get_name());
            return true;
        }
        bool rc = true;
        rc = evalate_extinct_taxa(otCLI, treeIndex, *raw, extinct_tips, relExtinctIDs, relIncSedIDs);
        return rc;
    }


    bool evalate_extinct_taxa(OTCLI & otCLI,
                             std::size_t treeIndex,
                             const TreeMappedWithSplits & tree,
                             const ConstNdPtrSet  & extinct_tips,
                             const OttIdSet & relevantExtinctIDs,
                             const OttIdSet & relevantIncSedIDs) {
        RTreeOttIDMapping<RTSplits> & taxdata = taxonomy->get_data();
        auto phylo_root = tree.get_root();
        const RTSplits & phylo_root_data = phylo_root->get_data();
        const auto & phylo_tip_ids = phylo_root_data.des_ids;
        OttIdSet fossil_ids;
        ConstNdPtrSet nonExtinctPhyloTips;
        OttIdSet non_fossil_ids;
        db_write_ott_id_set("phylo root des_ids", phylo_tip_ids);
        Phylo2Taxo phylo2taxo;
        _fill_initial_mapping(tree, extinct_tips, phylo_tip_ids, fossil_ids, nonExtinctPhyloTips, non_fossil_ids, phylo2taxo);
        ConstNdPtrSet taxoTips;
        for (auto i : phylo2taxo) {
            taxoTips.insert(i.second);
        }
        auto ind_tree_and_root = gen_taxo_induced_tree(phylo_tip_ids, phylo2taxo);
        const TaxoNd2IndParIDset & taxo_induced_tree = ind_tree_and_root.first;
        ConstNdPtr taxo_root = ind_tree_and_root.second;
        const std::size_t num_extinct_tips = extinct_tips.size();
        const std::size_t num_extant_tips = nonExtinctPhyloTips.size();
        const std::size_t total_num_tips = num_extinct_tips + num_extant_tips;
        if (taxo_induced_tree.size() == total_num_tips) {
            LOG(DEBUG) << "taxo_induced_tree.size() == total_num_tips EARLY EXIT\n";
            // if all if there is no structure in the induced tree, there cannot be any taxa contested..
            return true;
        }
        ConstNdPtrSet nonExtinctTaxoTips;
        for (auto extantPhyloTip : nonExtinctPhyloTips) {
            nonExtinctTaxoTips.insert(phylo2taxo.at(extantPhyloTip));
        }
        TaxoNd2RestrictedIDSet taxo_internal2id_set;
        for (auto taxoIndNdIt : taxo_induced_tree) {
            if (!contains(taxo_internal2id_set, taxoIndNdIt.second.first)) {
                taxo_internal2id_set[taxoIndNdIt.second.first] = taxoIndNdIt.second.second; 
            }
        }
        ConstNdPtrSet contestedByExtant;
        ConstNdPtrSet uncontestedByExtant;
        const OttIdSet extant_ids = set_difference_as_set(phylo_tip_ids, relevantExtinctIDs);
        for (auto tax2id_set_it : taxo_internal2id_set) {
            if (taxon_extant_conflicts_with_tree(tree, tax2id_set_it.second, extant_ids, relevantIncSedIDs)) {
                contestedByExtant.insert(tax2id_set_it.first);
            } else {
                uncontestedByExtant.insert(tax2id_set_it.first);
            }
        }
        std::string tree_id = tree.get_name();
        for (auto extinct_tip : extinct_tips) {
            TreeIDToPlacement & treeIDtoPlacement = unjoinedExtinctTaxonToPlacementSummary[extinct_tip->get_ott_id()];
            auto tree_name = tree.get_name();
            treeIDtoPlacement[tree_name] = evaluate_extinct_leaf_placement(tree, extinct_tip, extant_ids,
                                                                           relevantIncSedIDs, phylo2taxo, taxo_induced_tree,
                                                                           contestedByExtant, taxo_root);
        }
        return true;
    }

    NeedsMoveTipmostTaxonPair evaluate_extinct_leaf_placement(const TreeMappedWithSplits & tree,
                                                              ConstNdPtr extinct_tip,
                                                              const OttIdSet & extant_ids,
                                                              const OttIdSet & relevantIncSedIDs,
                                                              const Phylo2Taxo & phylo2taxo,
                                                              const TaxoNd2IndParIDset taxo_induced_tree,
                                                              const ConstNdPtrSet & contestedByExtant,
                                                              const ConstNdPtr taxo_root) {
        // A bit tricky here... as we move from the extinct tip rootward in each tree, we know
        //    that the one extinct tip we are analyzing is in each node's des_ids. So we don't bother
        //    adding it.
        NeedsMoveTipmostTaxonPair nmttp = {false, nullptr};
        OttId curr_id = extinct_tip->get_ott_id();
        ConstNdPtr rel_forking_phylo_anc = extinct_tip;
        OttIdSet rfpa_des = set_intersection_as_set(extant_ids, rel_forking_phylo_anc->get_data().des_ids);
        while (rfpa_des.size() < 1) {
            rel_forking_phylo_anc = rel_forking_phylo_anc->get_parent();
            if (rel_forking_phylo_anc == nullptr) {
                return nmttp;
            }
            const auto & np_des_ids = rel_forking_phylo_anc->get_data().des_ids;
            rfpa_des = set_intersection_as_set(extant_ids, np_des_ids);
        }
        if (rfpa_des == extant_ids) {
            return nmttp;
        }
        auto taxo_node = phylo2taxo.at(extinct_tip);
        const InducedParAndIds * tax_par_and_id_set = &(taxo_induced_tree.at(taxo_node));
        OttIdSet tax_forking_anc_des_ids;
        ConstNdPtr curr_taxo_nd = tax_par_and_id_set->first;
        while (true) {
            while (contains(contestedByExtant, curr_taxo_nd)){
                tax_par_and_id_set = &(taxo_induced_tree.at(curr_taxo_nd));
                curr_taxo_nd = tax_par_and_id_set->first;
            }
            if (curr_taxo_nd == taxo_root) {
                break;
            }
            tax_forking_anc_des_ids = set_intersection_as_set(tax_par_and_id_set->second, extant_ids);
            if (tax_forking_anc_des_ids.size() >= rfpa_des.size()) {
                if (is_subset(rfpa_des, tax_forking_anc_des_ids)) {
                    break;
                } else {
                    nmttp.first = true;
                    tax_par_and_id_set = &(taxo_induced_tree.at(curr_taxo_nd));
                    curr_taxo_nd = tax_par_and_id_set->first;
                }
            } else {
                if (!is_subset(tax_forking_anc_des_ids, rfpa_des)) {
                    assert(false); // MTH thinks this only happens if the extant taxa contest...
                    nmttp.first = true;
                }
                tax_par_and_id_set = &(taxo_induced_tree.at(curr_taxo_nd));
                curr_taxo_nd = tax_par_and_id_set->first;
            }
        }
        if (nmttp.first) {
            nmttp.second = curr_taxo_nd;
        }
        return nmttp;
    }

    bool taxon_extant_conflicts_with_tree(const TreeMappedWithSplits & tree,
                                          const OttIdSet & include_set,
                                          const OttIdSet & extant_ids,
                                          const OttIdSet & inc_sed_ids) {
        db_write_ott_id_set("taxon_extant_conflicts_with_tree", include_set);
        const OttIdSet extant_include = set_intersection_as_set(include_set, extant_ids);
        if (extant_include.size() < 2) {
            LOG(DEBUG) << "no/trivial include set....skipping";
        }
        const OttIdSet extant_exclude = set_difference_as_set(extant_ids, extant_include);
        const OttIdSet filtered_exclude = set_difference_as_set(extant_exclude, inc_sed_ids);
        if (filtered_exclude.empty()) {
            LOG(DEBUG) << "no exclude set....skipping";
        }
        for (auto phylo_internal : iter_node_internal_const(tree)) {
            const auto & phylo_inc_set = phylo_internal->get_data().des_ids;
            if (are_disjoint(phylo_inc_set, extant_include) || is_subset(extant_include, phylo_inc_set)) {
                continue;
            }
            if (have_intersection(phylo_inc_set, filtered_exclude)) {
                return true;
            }
        }
        return false;
    }
    
    using InducedTreeAndRoot = std::pair<TaxoNd2IndParIDset, ConstNdPtr>;
    InducedTreeAndRoot gen_taxo_induced_tree(const OttIdSet & phylo_tip_ids, const Phylo2Taxo & phylo2taxo) {
        TaxoNd2IndParIDset to_induced_mapping;
        ConstNdPtr root = nullptr;
        for (auto pt : phylo2taxo) {
            LOG(DEBUG) << "currPhylo = " << pt.first->get_name() << "\n";
            auto currTaxo = pt.second;
            assert(currTaxo != nullptr);
            OttIdSet currTaxoDes = set_intersection_as_set(phylo_tip_ids, currTaxo->get_data().des_ids);
            auto nextTaxo = currTaxo;
            OttIdSet nextTaxoDes;
            while (true) {
                nextTaxo = nextTaxo->get_parent();
                LOG(DEBUG) << "curr = " << currTaxo->get_name() << " next = " << nextTaxo->get_name();
                assert(nextTaxo != nullptr);
                nextTaxoDes = set_intersection_as_set(phylo_tip_ids, nextTaxo->get_data().des_ids);
                if (nextTaxoDes == currTaxoDes) {
                    LOG(DEBUG) << "...non-forking on induced";
                    continue;
                }
                LOG(DEBUG) << "...storing " << currTaxo->get_name() << " -> " << nextTaxo->get_name();
                to_induced_mapping[currTaxo] = InducedParAndIds(nextTaxo, nextTaxoDes);
                if (nextTaxoDes == phylo_tip_ids) {
                    assert(root == nullptr || root == nextTaxo);
                    root = nextTaxo;
                    LOG(DEBUG) << "...at root";
                    break;
                }
                currTaxo = nextTaxo;
                if (contains(to_induced_mapping, currTaxo)) {
                    LOG(DEBUG) << "...Already added";
                    break;
                }
                currTaxoDes = nextTaxoDes;
            }
        }
        assert(root != nullptr);
        return InducedTreeAndRoot(to_induced_mapping, root);
    }
    
    // fills fossil_ids, nonextinct_tips, non_fossil_ids, and phylo2taxo
    void _fill_initial_mapping(const TreeMappedWithSplits & tree,
                               const ConstNdPtrSet  & extinct_tips,
                               const OttIdSet & phylo_tip_ids,
                               OttIdSet & fossil_ids, 
                               ConstNdPtrSet & nonextinct_tips,
                               OttIdSet & non_fossil_ids, 
                               Phylo2Taxo & phylo2taxo
                               ) {
        RTreeOttIDMapping<RTSplits> & taxdata = taxonomy->get_data();
        for (auto nd : iter_leaf_const(tree)) {
            auto ottId = nd->get_ott_id();
            auto & nd_des_ids = nd->get_data().des_ids;
            auto tax_nd_for_tip = taxdata.ott_id_to_node.at(ottId);
            phylo2taxo[nd] = tax_nd_for_tip;
            auto & tax_nd_des_ids = tax_nd_for_tip->get_data().des_ids;
            auto tax_tip_des_overlap = set_intersection_as_set(tax_nd_des_ids, phylo_tip_ids);
            if (tax_tip_des_overlap.size() != 1) {
                db_write_ott_id_set("evidence of overlapping Ids: ", tax_tip_des_overlap);
                LOG(ERROR) << "Node \"" << nd->get_name() << "\" has taxonomic descendants that overlap with other tips in the tree.\n";
                throw OTCError("Non disjunct taxa as tips of tree");
            }
            if (contains(extinct_tips, nd)) {
                fossil_ids.insert(ottId);
            } else {
                non_fossil_ids.insert(ottId);
                nonextinct_tips.insert(nd);
            }
        }
    }
    
};


bool handleJSONOutput(OTCLI & otCLI, const std::string &narg);
bool handleExtinctJSON(OTCLI & , const std::string &);
bool handleExtinctJSONInline(OTCLI & , const std::string &) ;


bool handleJSONOutput(OTCLI & otCLI, const std::string & narg) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting filepath of output JSON after the -j argument.");
    }
    proc->outputJSONLogFile = narg;
    proc->useStdOut = false;
    return true;
}


bool handleExtinctJSON(OTCLI & otCLI, const std::string & narg) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting filepath to the exinction ID JSON file  the -t argument.");
    }
    proc->extinctInputJSONFilepath = narg;
    return true;
}

bool handleExtinctJSONInline(OTCLI & otCLI, const std::string & narg) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting JSON content inline after the -i argument.");
    }
    proc->extinctInputJSONInline = narg;
    proc->useCmdLine = true;
    return true;
}


int main(int argc, char *argv[]) {
    const char * helpMsg =  "This tool treats the fossil taxa in the input taxonomic tree as poorly placed. " \
        "It will produce a set of edits (encoded as JSON) that can be applied to the taxonomy. If these " \
        "edits were applied to the taxonomy, the none of the fossil taxa included in the phylogenetic trees " \
        "would increase the set of taxa that are contested by input trees.\n" \
        "Input is a full taxonomy tree some number of input trees, and a (JSON-formatted) list of taxon IDs that are extinct.\n" \
        "Extinct taxon info can be given as a file name (-t flag) or on the command-line (-i flag)\n" \
        "The output JSON can be specified as stdout (default) or as filename (-j flag)\n";
        
    OTCLI otCLI("otc-move-extinct-higher-to-avoid-contesting-taxa",
                helpMsg,
                "-jedits.json -textinct.json taxonomy.tre inp1.tre inp2.tre");
    MoveExtinctHigherState proc;
    otCLI.add_flag('j',
                  "Name of an output JSON file that summarizes the set of IDs used to exemplify each taxon.",
                  handleJSONOutput,
                  true);
    otCLI.add_flag('t',
                  "Name of input JSON file listing the exinct OTT IDs",
                  handleExtinctJSON,
                  true);
    otCLI.add_flag('i',
                  "exinct OTT IDs as JSON argument, not the filepath to JSON",
                  handleExtinctJSONInline,
                  true);
    return tax_dependent_tree_processing_main(otCLI, argc, argv, proc, 2, true);
}
