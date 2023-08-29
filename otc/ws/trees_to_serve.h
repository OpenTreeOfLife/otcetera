#ifndef TREES_TO_SERVE_H
#define TREES_TO_SERVE_H

#include <set>
#include <thread>
#include "otc/ws/tolws.h"
#include "otc/taxonomy/patching.h"

namespace otc
{
using ReadableTaxonomy = std::pair<const RichTaxonomy &,
                                   std::unique_ptr<ReadMutexWrapper> >;
using WritableTaxonomy = std::pair<RichTaxonomy &,
                                   std::unique_ptr<WriteMutexWrapper> >;
                                   
class TreesToServe {
    std::list< SummaryTreeAnnotation> annotation_list;
    std::list<std::unique_ptr<SummaryTree_t> > tree_list;
    std::map<std::string, const SummaryTree_t *> id_to_tree;
    std::map<std::string, const SummaryTreeAnnotation *> id_to_annotations;
    std::string default_synth_id;
    PatchableTaxonomy * taxonomy_ptr = nullptr;
    std::map<std::string, const std::string *> stored_strings;
    std::list<std::string> stored_strings_list;
    const RichTaxTree * taxonomy_tree = nullptr;
    vec_src_node_id_mapper src_node_id_storer;
    std::map<src_node_id, std::uint32_t> lookup_for_node_ids_while_registering_trees;
    bool finalized = false;
    mutable ParallelReadSerialWrite taxonomy_thread_safety;

public:
    explicit TreesToServe();

    const src_node_id & decode_study_node_id_index(std::uint32_t sni_ind) const;

    std::uint32_t get_source_node_id_index(src_node_id sni);

    const std::string * get_stored_string(const std::string & k);

    void set_taxonomy(PatchableTaxonomy &taxonomy);

    using ReadableTaxonomy = std::pair<const PatchableTaxonomy &, std::unique_ptr<ReadMutexWrapper> >;
    using WritableTaxonomy = std::pair<PatchableTaxonomy &, std::unique_ptr<WriteMutexWrapper> >;
    ReadableTaxonomy get_readable_taxonomy() const;

    WritableTaxonomy get_writable_taxonomy();

    void fill_ott_id_set(const std::bitset<32> & flags,
                         OttIdSet & ott_id_set,
                         OttIdSet & suppressed_from_tree);

    typedef std::pair<SummaryTree_t &, SummaryTreeAnnotation &> SumTreeInitPair;
    SumTreeInitPair get_new_tree_and_annotations(const std::string & configfilename,
                                                 const std::string & filename);

    void register_last_tree_and_annotations();

    void free_last_tree_and_annotations();

    std::string get_default_tree() const;

    std::set<std::string> get_available_trees() const;

    const SummaryTreeAnnotation * get_annotations(std::string synth_id) const;

    const SummaryTree_t * get_summary_tree(std::string synth_id) const;

    std::size_t get_num_trees() const;

    void final_tree_added();
};



} // namespace otc
#endif
