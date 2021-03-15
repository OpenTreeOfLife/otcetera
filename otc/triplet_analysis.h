#ifndef OTCETERA_TRIPLET_ANALYSIS_H
#define OTCETERA_TRIPLET_ANALYSIS_H
#include <ostream>
#include <vector>
#include "otc/tree_as_split_set.h"
#include "otc/triple_dist.h"
#include "otc/all_triplets.h"

namespace otc {

template<typename T>
class TripletDistAnalysis {
    public:
    using node_t = typename T::node_type;
    using paired_nodes = std::pair<const node_t *, const node_t *>;
    using pair_indices = std::pair<std::size_t, std::size_t>;
    
    private:
    const T & tree1;
    const T & tree2;
    std::vector<pair_indices> tot_diff_comp_by_round;
    std::vector<paired_nodes> pruned_after_each_round;
    std::vector<pair_indices> highest_diff_comp_by_round;
    std::set<std::size_t> pruned_inds;

    public:
    TripletDistAnalysis(const T & inp_tre1,
                        const T & inp_tre2)
        :tree1(inp_tre1),
        tree2(inp_tre2) {
        this->run();
    }

    std::size_t get_num_rounds() const {
        return tot_diff_comp_by_round.size();
    }

    pair_indices get_tot_diff_comp_for_round(std::size_t ind) const {
        return tot_diff_comp_by_round.at(ind);
    }

    paired_nodes get_nodes_paired_after_round(std::size_t ind) const {
        if (ind >= highest_diff_comp_by_round.size()) {
            return paired_nodes{nullptr, nullptr};
        }
        return pruned_after_each_round.at(ind);
    }

    pair_indices get_highest_diff_comp_for_round(std::size_t ind) const {
        if (ind >= highest_diff_comp_by_round.size()) {
            return pair_indices{0, 0};
        }
        return highest_diff_comp_by_round.at(ind);
    }

    private:
    void run() {
        using TreeAsUIntSplits = GenTreeAsUIntSplits<T>;
        using AllConflictTriplets = AllTriplets<T>;
        TreeAsUIntSplits tas_1{tree1};
        TreeAsUIntSplits tas_2{tree2};
        if (tas_1.leaf_label_to_ind != tas_2.leaf_label_to_ind) {
            throw OTCError() << "trees must have the same leaf label set.\n";
        }
        AllConflictTriplets t_1_rt{tas_1};
        AllConflictTriplets t_2_rt{tas_2};
        TripletDist rtdist{t_1_rt, t_2_rt};
        std::size_t most_rec_pruned = 0;
        for (;;) {
            const auto dc = rtdist.calc_diff_comp(pruned_inds);
            tot_diff_comp_by_round.push_back(dc);
            if (dc.first < 1) {
                break;
            }
            const auto hds = rtdist.get_highest_dist(pruned_inds);
            assert(!hds.empty());
            most_rec_pruned = *(hds.begin());
            pruned_inds.insert(most_rec_pruned);
            auto ndp1 = tas_1.ind_to_nd.at(most_rec_pruned);
            auto ndp2 = tas_2.ind_to_nd.at(most_rec_pruned);
            pruned_after_each_round.push_back(paired_nodes{ndp1, ndp2});
            auto hdc = rtdist.get_diff_comp_for_index(most_rec_pruned);
            highest_diff_comp_by_round.push_back(hdc);
        }
    }
};


} // namespace

#endif
