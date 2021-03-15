#ifndef OTCETERA_TRIPLE_DIST_H
#define OTCETERA_TRIPLE_DIST_H
#include <ostream>
#include <vector>
#include "otc/otc_base_includes.h"
#include "otc/all_triplets.h"
#include "otc/util.h"

namespace otc {


enum RT_COMP {BOTH_UNRES = 0,
              COMPAT = 1,
              SAME_RES = 2,
              CONFLICT_RES = 3,
              NO_COMP = 4
            };

RT_COMP comp_tt(const TRIPLET_TYPE & qt1, const TRIPLET_TYPE & qt2);

inline RT_COMP comp_tt(const TRIPLET_TYPE & qt1, const TRIPLET_TYPE & qt2) {
    if (qt1 == TRIPLET_TYPE::UNKNOWN || qt2 == TRIPLET_TYPE::UNKNOWN
        || qt1 == TRIPLET_TYPE::NOT_Q || qt2 == TRIPLET_TYPE::NOT_Q) {
        return RT_COMP::NO_COMP;
    }
    if (qt1 == TRIPLET_TYPE::POLYTOMY) {
        if (qt2 == TRIPLET_TYPE::POLYTOMY) {
            return RT_COMP::BOTH_UNRES;
        }
        return RT_COMP::COMPAT;
    }
    if (qt1 == qt2) {
        return RT_COMP::SAME_RES;
    }
    return RT_COMP::CONFLICT_RES;
}


template <typename T>
class TripletDist {
    const AllTriplets<T> & rttree1;
    const AllTriplets<T> & rttree2;
    std::size_t num_tips;
    std::size_t round;
    std::size_t num_diffs;
    std::size_t num_comp;

    std::vector<std::size_t> diff_by_taxon;
    std::vector<std::size_t> comp_by_taxon;
    
    using off_by_2_cmp_v = std::vector<RT_COMP> ;
    using off_by_1_cpm_v = std::vector<off_by_2_cmp_v>;
    using top_level_cmp = std::vector<off_by_1_cpm_v>;
    
    top_level_cmp by_lowest;
    std::set<std::size_t> prev_ignored;
    public:
    using pair_indices = std::pair<std::size_t, std::size_t>;

    TripletDist(const AllTriplets<T> & q1,
                const AllTriplets<T> & q2)
    :rttree1(q1),
     rttree2(q2),
     round(0) {
        num_tips = q1.num_tips;
        assert(num_tips == q2.num_tips);
        std::set<std::size_t> empty_set;
        this->calc_diffs_mat(empty_set, true);
    }

    pair_indices get_diff_comp() const {
        return {num_diffs, num_comp};
    }

    pair_indices get_diff_comp_for_index(std::size_t ind) const {
        return {diff_by_taxon.at(ind), comp_by_taxon.at(ind)};
    }


    std::set<std::size_t> get_highest_dist(std::set<std::size_t> & taboo) const {
        std::set<std::size_t> biggest_dist_inds;
        double frac_diff = -1;
        for (std::size_t i = 0; i < diff_by_taxon.size(); ++i) {
            if (contains(taboo, i)) {
                continue;
            }
            const pair_indices num_diff{diff_by_taxon.at(i), comp_by_taxon.at(i)};
            const auto cd = frac_diff_from_pair(num_diff);
            if (cd > frac_diff) {
                biggest_dist_inds.clear();
                biggest_dist_inds.insert(i);
                frac_diff = cd;
            } else if (cd == frac_diff) {
                biggest_dist_inds.insert(i);
            }
        }
        return biggest_dist_inds;
    }


    pair_indices calc_diff_comp(const std::set<std::size_t> & taboo) {
        if (taboo != prev_ignored) {
            this->calc_diffs_mat(taboo, false);
        }
        return this->get_diff_comp();
    }

    private:
    void calc_diffs_mat(const std::set<std::size_t> & taboo, 
                        bool calc_comp) {
        prev_ignored = taboo;
        std::size_t n_r = static_cast<std::size_t>(static_cast<int>(num_tips) - 2);
        by_lowest.reserve(n_r);
        const RT_COMP def = RT_COMP::NO_COMP;
        for (std::size_t row_n = 0; row_n < n_r; ++row_n) {
            by_lowest.push_back(gen_tm_by_sec<RT_COMP>(num_tips, row_n, def));
        }
        diff_by_taxon.assign(num_tips, 0);
        comp_by_taxon.assign(num_tips, 0);
        const auto & bl1 = rttree1.by_lowest;
        const auto & bl2 = rttree2.by_lowest;
        const std::size_t ie = num_tips - 2;
        const std::size_t je = num_tips - 1;
        const std::size_t ke = num_tips;
        num_diffs = num_comp = 0;
        for (std::size_t i = 0; i < ie; ++i) {
            if (contains(taboo, i)) {
                continue;
            }
            const auto & ibl1 = bl1.at(i);
            const auto & ibl2 = bl2.at(i);
            auto & idbl = by_lowest.at(i);
            for (std::size_t j = i + 1; j < je; ++j) {
                if (contains(taboo, j)) {
                    continue;
                }
                const auto rel_j = j - i - 1; 
                const auto & jbl1 = ibl1.at(rel_j);
                const auto & jbl2 = ibl2.at(rel_j);
                auto & jdbl = idbl.at(rel_j);
                for (std::size_t k = j + 1; k < ke; ++k) {
                    if (contains(taboo, k)) {
                        continue;
                    }
                    const auto rel_k = k - j - 1;
                    const auto & el1 = jbl1.at(rel_k);
                    const auto & el2 = jbl2.at(rel_k);
                    RT_COMP qcmp = RT_COMP::NO_COMP;
                    if (calc_comp) {
                        qcmp = comp_tt(el1, el2);
                        jdbl.at(rel_k) = qcmp;
                    } else {
                        qcmp = jdbl.at(rel_k);
                    }
                    if (qcmp == RT_COMP::NO_COMP) {
                        continue;
                    }
                    // here we'll just count conflicts as distances,
                    //  so polytomy, compat or same all count as no diff
                    if (qcmp == RT_COMP::CONFLICT_RES) {
                        diff_by_taxon.at(i) += 1;
                        diff_by_taxon.at(j) += 1;
                        diff_by_taxon.at(k) += 1;
                        num_diffs += 1;
                    }
                    num_comp += 1;
                    comp_by_taxon.at(i) += 1;
                    comp_by_taxon.at(j) += 1;
                    comp_by_taxon.at(k) += 1;
                }
            }
        }
    }
};



} // namespace otc

#endif
