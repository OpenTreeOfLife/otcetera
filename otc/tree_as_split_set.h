#ifndef OTCETERA_TREES_AS_SPLIT_SETS_H
#define OTCETERA_TREES_AS_SPLIT_SETS_H
#include <map>
#include <string>
#include <vector>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/tree.h"
#include "otc/tree_iter.h"
#include "otc/util.h"
namespace otc {

template<typename T>
class GenTreeAsSplits {
    public:
    using node_t = typename T::node_type;
    using str_set = std::set<std::string>;

    std::map<std::string, const node_t *> leaf_label_to_nd;
    std::map<const node_t *, str_set> nd_to_taxset;
    std::map<str_set, const node_t *> inf_taxset_to_nd;
    str_set leaf_labels;

    GenTreeAsSplits(const T & tree) {
        auto nodes = all_nodes_postorder(tree);
        for (auto nd : nodes) {
            str_set ts;
            if (nd->is_tip()) {
                auto & nl = nd->get_name();
                ts.insert(nl);
                leaf_labels.insert(nl);
                leaf_label_to_nd[nl] = nd;
            } else {
                for (auto c : iter_child_const(*nd)) {
                    auto cts = nd_to_taxset[c];
                    ts.insert(cts.begin(), cts.end());
                }
                if (nd->get_parent() != nullptr) {
                    inf_taxset_to_nd[ts] = nd;
                }
            }
            nd_to_taxset[nd] = ts;
        }
    }
};

template<typename T>
class GenTreeAsUIntSplits {
    public:
    using node_t = typename T::node_type;
    using uint_set = std::set<std::size_t>;

    std::map<std::string, std::size_t> leaf_label_to_ind;
    std::vector<const node_t *> ind_to_nd;
    std::map<const node_t *, uint_set> nd_to_taxset;
    std::map<uint_set, const node_t *> inf_taxset_to_nd;

    GenTreeAsUIntSplits(const T & tree) {
        const GenTreeAsSplits<T> str_tas{tree};
        std::size_t ind = 0;
        this->ind_to_nd.reserve(str_tas.leaf_labels.size());
        for (auto label : str_tas.leaf_labels) {
            this->leaf_label_to_ind[label] = ind;
            auto nd_ptr = str_tas.leaf_label_to_nd.at(label);
            this->ind_to_nd.push_back(nd_ptr);
            ++ind;
        }
        for (auto labs_nd_it : str_tas.inf_taxset_to_nd) {
            const auto & label_set = labs_nd_it.first;
            const node_t * nd_ptr = labs_nd_it.second;
            uint_set ind_set;
            for (auto label: label_set) {
                auto lab_ind = this->leaf_label_to_ind.at(label);
                ind_set.insert(lab_ind);
            }
            this->nd_to_taxset[nd_ptr] = ind_set;
            this->inf_taxset_to_nd[ind_set] = nd_ptr;
        }
    }
};



} // namespace otc
#endif
