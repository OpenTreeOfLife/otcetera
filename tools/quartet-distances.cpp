// See https://github.com/OpenTreeOfLife/opentree/wiki/Open-Tree-of-Life-APIs-v3#synthetic-tree
#include <iostream>
#include <tuple>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include "json.hpp"
#include "otc/conflict.h"
#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/tree_as_split_set.h"


using json=nlohmann::json;


using namespace otc;
using std::vector;
template <typename T, typename U>
using Map = std::unordered_multimap<T,U>;
template <typename T, typename U>
using map = std::unordered_map<T,U>;
template <typename T>
using set = std::unordered_set<T>;
using std::string;
//using std::map;
using std::pair;
using std::tuple;
using std::unique_ptr;
using std::string;

namespace po = boost::program_options;
using po::variables_map;

// TODO: Could we exemplify tips here, if we had access to the taxonomy?
variables_map parse_cmd_line(int argc,char* argv[]) 
{ 
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("first", value<string>(),"First tree filepath.")
        ("second", value<string>(),"Second tree filepath.")
        ;

    options_description reporting("Reporting options");
    // reporting.add_options()
    //     ("each",value<bool>()->default_value(false),"Show separate results for each input tree")
    //     ("all",value<bool>()->default_value(true),"Show accumulated over all input trees")
    //     ("switch","Count synth nodes instead of input tree nodes")
    //     ("names,N","Write out node names instead of counts.")
    //     ;

    options_description visible;
    visible.add(reporting).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("first", 1);
    p.add("second", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-quartet-distances <first> <second> [OPTIONS]\n"
                                                    "Reports stats based on tree differences using quartet distances.",
                                                    visible, invisible, p);

    return vm;
}

namespace std {

template<typename T, typename U>
struct hash<pair<T,U>> {
    std::size_t operator()(const std::pair<T,U>& p) const noexcept {return std::hash<T>()(p.first) * std::hash<U>()(p.second);}
};

} // namespace std

using Tree_t = ConflictTree;
using node_t = Tree_t::node_type;

using str_set = std::set<std::string>;
using TreeAsUIntSplits = GenTreeAsUIntSplits<Tree_t>;

enum QUARTET_TYPE {UNKNOWN = 0,
                   POLYTOMY = 1,
                   ONE_TWO = 2,
                   ONE_THREE = 3,
                   ONE_FOUR = 4,
                   NOT_Q = 5
                   };

enum Q_COMP {BOTH_UNRES = 0,
             COMPAT = 1,
             SAME_RES = 2,
             CONFLICT_RES = 3,
             NO_COMP = 4
            };

inline Q_COMP comp_qt(const QUARTET_TYPE & qt1, const QUARTET_TYPE & qt2) {
    if (qt1 == QUARTET_TYPE::UNKNOWN || qt2 == QUARTET_TYPE::UNKNOWN
        || qt1 == QUARTET_TYPE::NOT_Q || qt2 == QUARTET_TYPE::NOT_Q) {
        return Q_COMP::NO_COMP;
    }
    if (qt1 == QUARTET_TYPE::POLYTOMY) {
        if (qt2 == QUARTET_TYPE::POLYTOMY) {
            return Q_COMP::BOTH_UNRES;
        }
        return Q_COMP::COMPAT;
    }
    if (qt1 == qt2) {
        return Q_COMP::SAME_RES;
    }
    return Q_COMP::CONFLICT_RES;
}

inline void write_qt(std::ostream & out, const QUARTET_TYPE & qt) {
    switch(qt) {
        case QUARTET_TYPE::UNKNOWN : out << "?" ; break;
        case QUARTET_TYPE::POLYTOMY : out << "*" ; break;
        case QUARTET_TYPE::NOT_Q : out << "X" ; break;
        case QUARTET_TYPE::ONE_TWO : out << "12" ; break;
        case QUARTET_TYPE::ONE_THREE : out << "13" ; break;
        case QUARTET_TYPE::ONE_FOUR : out << "14" ; break;
    }
}


template<typename T>
inline std::vector<T> 
gen_by_fourth(std::size_t num_tax, std::size_t third_index, const T & def) {
    //std::cerr << "    gen_by_fourth(" << num_tax << ", " << third_index << ")\n";
    const std::size_t min_real_ind = third_index + 1;
    assert(min_real_ind < num_tax);
    // N - 1 (for the zero-based ind)
    const std::size_t max_real_ind = num_tax - 1;
    //  std::cerr << "  max_real_ind=" << max_real_ind << " min_real_ind=" << min_real_ind << "\n";
    assert(max_real_ind >= min_real_ind);
    const std::size_t my_size = 1 + max_real_ind - min_real_ind;
    std::vector<T> ret{my_size, def};
    return ret;
}

template<typename T>
inline std::vector<std::vector<T> >
gen_by_third(std::size_t num_tax, std::size_t sec_index, const T & def ) {
    //std::cerr << "  gen_by_third(" << num_tax << ", " << sec_index << ")\n";
    std::vector<std::vector<T> > ret;
    const std::size_t min_real_ind = sec_index + 1;
    assert(min_real_ind < num_tax);
    // N - 1 (for 2 more indices) - 1 (for the zero-based ind)
    const std::size_t max_real_ind = num_tax - 2;
    // std::cerr << "  max_real_ind=" << max_real_ind << " min_real_ind=" << min_real_ind << "\n";
    assert(max_real_ind >= min_real_ind);
    const std::size_t my_size = 1 + max_real_ind - min_real_ind;
    ret.reserve(my_size);
    for (std::size_t my_real_ind = min_real_ind; my_real_ind <=  max_real_ind; ++my_real_ind) {
        ret.push_back(gen_by_fourth(num_tax, my_real_ind, def));
    }
    return ret;
}

template<typename T>
inline std::vector<std::vector<std::vector<T> > >
gen_by_sec(std::size_t num_tax, std::size_t first_ind, const T & def) {
    //std::cerr << "gen_by_sec(" << num_tax << ", " << first_ind << ")\n";
    std::vector<std::vector<std::vector<T> > > ret;
    const std::size_t min_real_ind = first_ind + 1;
    assert(min_real_ind < num_tax);
    // N - 2 (for 2 more indices) - 1 (for the zero-based ind)
    const std::size_t max_real_ind = num_tax - 3;
    // std::cerr << "  max_real_ind=" << max_real_ind << " min_real_ind=" << min_real_ind << "\n";
    assert(max_real_ind >= min_real_ind);
    const std::size_t my_size = 1 + max_real_ind - min_real_ind;
    ret.reserve(my_size);
    for (std::size_t my_real_ind = min_real_ind; my_real_ind <=  max_real_ind; ++my_real_ind) {
        ret.push_back(gen_by_third(num_tax, my_real_ind, def));
    }
    return ret;
}

class AllQuartets {
    using off_by_3_v = std::vector<QUARTET_TYPE> ;
    using off_by_2_v = std::vector<off_by_3_v>;
    using off_by_1_v = std::vector<off_by_2_v>;
    using top_level = std::vector<off_by_1_v>;

    top_level by_lowest;

    public:
    using uint_set = std::set<std::size_t>;
    std::size_t num_tips;

    AllQuartets(const TreeAsUIntSplits & tas) {
        num_tips = tas.ind_to_nd.size();
        if (num_tips < 4) {
            return;
        }
        std::size_t n_r = static_cast<std::size_t>(static_cast<int>(num_tips) - 3);
        by_lowest.reserve(n_r);
        const QUARTET_TYPE def = QUARTET_TYPE::UNKNOWN;
        for (std::size_t row_n = 0; row_n < n_r; ++row_n) {
            by_lowest.push_back(gen_by_sec<QUARTET_TYPE>(num_tips, row_n, def));
        }
        
        uint_set full_ind_set;
        for (std::size_t i = 0; i < num_tips; ++i) {
            full_ind_set.insert(i);
        }
        for (auto inf_indset_nd_pair : tas.inf_taxset_to_nd) {
            const uint_set & inf_ind_set = inf_indset_nd_pair.first;
            const auto outgroup = set_difference_as_set(full_ind_set, inf_ind_set);
            const node_t * par_nd = inf_indset_nd_pair.second;
            this->register_nd(par_nd, outgroup, tas);
        }
        uint_set empty_set;
        this->register_nd(tas.root, empty_set, tas);
    }


    void write(std::ostream & out) const {
        std::size_t i = 0;
        // out << by_lowest.size() << " top_level rows\n";
        for (const auto & row : by_lowest) {
            write_by_sec(out, row, i++);
        }
    }
    private:

    void register_nd(const node_t * par_nd, 
                     const uint_set & outgroup,
                     const TreeAsUIntSplits & tas) {
        const node_t * curr_ch = par_nd->get_first_child();
        while (curr_ch != nullptr) {
            const auto & cc_ind_set = tas.nd_to_taxset.at(curr_ch);
            const node_t * curr_sib = curr_ch->get_next_sib();
            while (curr_sib != nullptr) {
                const auto & cs_ind_set = tas.nd_to_taxset.at(curr_sib);
                this->register_sibs(cc_ind_set, cs_ind_set, outgroup);
                curr_sib = curr_sib->get_next_sib();
            }
            curr_ch = curr_ch->get_next_sib();
        }
        const auto nd_out_deg = par_nd->get_out_degree();
        if (nd_out_deg > 2) {
            register_polytomy(par_nd, outgroup, tas);
        }
    }

    void register_polytomy(const node_t * par_nd, 
                           const uint_set & outgroup,
                           const TreeAsUIntSplits & tas) {
        const auto & nd_to_taxset = tas.nd_to_taxset;
        const node_t * curr_ch = par_nd->get_first_child();
        while (curr_ch != nullptr) {
            const node_t * curr_f_sib = curr_ch->get_next_sib();
            const auto & cc_ind_set = nd_to_taxset.at(curr_ch);
            while (curr_f_sib != nullptr) {
                const auto & cf_ind_set = nd_to_taxset.at(curr_f_sib);
                const node_t * curr_s_sib = curr_f_sib->get_next_sib();
                while (curr_s_sib != nullptr) {
                    const auto & cs_ind_set = nd_to_taxset.at(curr_s_sib);
                    this->register_poly_out(cc_ind_set, cf_ind_set, cs_ind_set, outgroup);
                    const node_t * curr_t_sib = curr_s_sib->get_next_sib();
                    while (curr_t_sib != nullptr) {
                        const auto & ct_ind_set = nd_to_taxset.at(curr_t_sib);
                        this->register_poly_out(cc_ind_set, cf_ind_set, cs_ind_set, ct_ind_set);
                        curr_t_sib = curr_t_sib->get_next_sib();
                    }
                    curr_s_sib = curr_s_sib->get_next_sib();
                }
                curr_f_sib = curr_f_sib->get_next_sib();
            }
            curr_ch = curr_ch->get_next_sib();
        }
    }
    void register_poly_out(const uint_set & f_ind_set, 
                       const uint_set & s_ind_set,
                       const uint_set & t_ind_set,
                       const uint_set & out_ind_set) {
        std::size_t fs_small, fs_large;
        std::size_t fst_small, fst_mid, fst_large;
        for (auto fci : f_ind_set) {
            for (auto sci : s_ind_set) {
                fs_small = fci;
                fs_large = sci;
                if (fs_small > fs_large) {
                    std::swap(fs_small, fs_large);
                }
                for (auto tci : t_ind_set) {
                    if (tci < fs_small) {
                        fst_small = tci;
                        fst_mid = fs_small;
                        fst_large = fs_large;
                    } else if (tci < fs_large) {
                        fst_small = fs_small;
                        fst_mid = tci;
                        fst_large = fs_large;
                    } else {
                        fst_small = fs_small;
                        fst_mid = fs_large;
                        fst_large = tci;
                    }
                    for (auto oci : out_ind_set) {
                        this->register_poly_last_unsorted(fst_small, fst_mid, fst_large, oci);
                    }
                }
            }
        }   
    }

    void register_poly_last_unsorted(std::size_t u1, std::size_t u2,
                          std::size_t u3, std::size_t uu) {
        std::size_t s1, s2, s3, s4;
        if (uu < u2) {
            s3 = u2; s4 = u3;
            if (uu < u1) {
                s1 = uu; s2 = u1; 
            } else {
                s1 = u1; s2 = uu; 
            }
        } else {
            s1 = u1; s2 = u2;
            if (uu < u3) {
                s3 = uu ; s4 = u3;
            } else {
                s3 = u3 ; s4 = uu;
            }
        }
        this->register_sorted(QUARTET_TYPE::POLYTOMY, s1, s2, s3, s4);
    }

    void register_sibs(const uint_set & lc_ind_set, 
                       const uint_set & nc_ind_set,
                       const uint_set & out_ind_set) {
        std::size_t in_small, in_large;
        for (auto lci : lc_ind_set) {
            for (auto nci : nc_ind_set) {
                if (lci < nci) {
                    in_small = lci;
                    in_large = nci;
                } else {
                    in_small = nci;
                    in_large = lci;
                }
                for (auto oi_it = out_ind_set.begin(); oi_it != out_ind_set.end(); ++oi_it) {
                    auto noi_it = oi_it;
                    for (++noi_it; noi_it != out_ind_set.end(); ++noi_it) {
                        register_quartet(in_small, in_large, *oi_it, *noi_it);
                    }
                }
            }
        }
    }

    void register_quartet(std::size_t in_small, std::size_t in_large,
                          std::size_t out_small, std::size_t out_large) {
        assert(in_small < in_large);
        assert(out_small < out_large);
        if (in_small < out_small) {
            if (in_large < out_small) {
                this->register_sorted(QUARTET_TYPE::ONE_TWO, in_small, in_large, out_small, out_large);
            } else if (in_large < out_large) {
                assert(in_large != out_small);
                this->register_sorted(QUARTET_TYPE::ONE_THREE, in_small, out_small, in_large, out_large);
            } else {
                assert(out_large < in_large);
                this->register_sorted(QUARTET_TYPE::ONE_FOUR, in_small, out_small, out_large, in_large);
            }
        } else if (in_small < out_large) {
            assert(out_small < in_small);
            if (in_large < out_large) {
                this->register_sorted(QUARTET_TYPE::ONE_FOUR, out_small, in_small, in_large, out_large);
            } else {
                assert(out_large < in_large);
                this->register_sorted(QUARTET_TYPE::ONE_THREE, out_small, in_small, out_large, in_large);
            }
        } else {
            assert(out_large < in_small);
            this->register_sorted(QUARTET_TYPE::ONE_TWO, out_small, out_large, in_small, in_large);
        }
    }

    void register_sorted(QUARTET_TYPE qt, std::size_t fir_ind, std::size_t sec_ind,
                          std::size_t thi_ind, std::size_t fou_ind) {
        assert(fou_ind > thi_ind);
        std::size_t rel_fou = fou_ind - thi_ind - 1;
        assert(thi_ind > sec_ind);    
        std::size_t rel_thi = thi_ind - sec_ind - 1;
        assert(sec_ind > fir_ind);
        std::size_t rel_sec = sec_ind - fir_ind - 1;
        try {
            this->by_lowest.at(fir_ind).at(rel_sec).at(rel_thi).at(rel_fou) = qt;
        } catch (...) {
            std::cerr << "Error registering ";
            write_qt(std::cerr, qt);
            std::cerr << " at (" << fir_ind << ", " << sec_ind << ", " << thi_ind << ", " << fou_ind << ")" << std::endl;
            throw;
        }
    }
    

    void write_by_sec(std::ostream & out, const off_by_1_v & by_sec, std::size_t first) const {
        std::size_t i = first + 1;
        // out << by_sec.size() << " 2nd-level rows\n";
        for (const auto & row : by_sec) {
            write_by_third(out, row, first, i++);
        }
    }

    void write_by_third(std::ostream & out, const off_by_2_v & by_third,
                        std::size_t first, std::size_t sec) const {
        std::size_t i = sec + 1;
        // out << by_third.size() << " 3rd-level rows\n";
        for (const auto & row : by_third) {
            write_by_fourth(out, row, first, sec, i++);
        }
    }

    void write_by_fourth(std::ostream & out, const off_by_3_v & by_fourth,
                         std::size_t first, std::size_t sec, std::size_t third) const {
        std::size_t i = third + 1;
        for (const QUARTET_TYPE el : by_fourth) {
            if (i > third + 1) {
                out << " ";
            }
            out << "q(" << first << ", " << sec << ", " << third << ", " << i << ")=";
            write_qt(out, el);
            i++;
        }
        out << '\n';
    }

    friend class QuartDist;
};

inline double frac_diff_from_pair(const std::pair<std::size_t, std::size_t> & p) {
    const double n = static_cast<double>(p.first);
    const double d = static_cast<double>(p.second);
    return n/d;
}

class QuartDist {
    const AllQuartets & qtree1;
    const AllQuartets & qtree2;
    std::size_t num_tips;
    std::size_t round;
    std::size_t num_diffs;
    std::size_t num_comp;

    std::vector<std::size_t> diff_by_taxon;
    std::vector<std::size_t> comp_by_taxon;
    
    using off_by_3_cmp_v = std::vector<Q_COMP> ;
    using off_by_2_cmp_v = std::vector<off_by_3_cmp_v>;
    using off_by_1_cpm_v = std::vector<off_by_2_cmp_v>;
    using top_level_cmp = std::vector<off_by_1_cpm_v>;

    top_level_cmp by_lowest;
    public:

    QuartDist(const AllQuartets & q1,
              const AllQuartets & q2)
    :qtree1(q1),
    qtree2(q2),
    round(0) {
        num_tips = q1.num_tips;
        assert(num_tips == q2.num_tips);
        this->calc_diffs_mat();
    }

    std::pair<std::size_t, std::size_t> get_diff_comp() const {
        return {num_diffs, num_comp};
    }


    private:
    void calc_diffs_mat() {
        std::size_t n_r = static_cast<std::size_t>(static_cast<int>(num_tips) - 3);
        by_lowest.reserve(n_r);
        const Q_COMP def = Q_COMP::NO_COMP;
        for (std::size_t row_n = 0; row_n < n_r; ++row_n) {
            by_lowest.push_back(gen_by_sec<Q_COMP>(num_tips, row_n, def));
        }

        diff_by_taxon.assign(num_tips, 0);
        comp_by_taxon.assign(num_tips, 0);
        const auto & bl1 = qtree1.by_lowest;
        const auto & bl2 = qtree2.by_lowest;
        const std::size_t ie = num_tips - 3;
        const std::size_t je = num_tips - 2;
        const std::size_t ke = num_tips - 1;
        const std::size_t le = num_tips;
        num_diffs = num_comp = 0;
        for (std::size_t i = 0; i < ie; ++i) {
            const auto & ibl1 = bl1.at(i);
            const auto & ibl2 = bl2.at(i);
            auto & idbl = by_lowest.at(i);

            for (std::size_t j = i + 1; j < je; ++j) {
                const auto rel_j = j - i - 1; 
                const auto & jbl1 = ibl1.at(rel_j);
                const auto & jbl2 = ibl2.at(rel_j);
                auto & jdbl = idbl.at(rel_j);

                for (std::size_t k = j + 1; k < ke; ++k) {
                    const auto rel_k = k - j - 1;
                    const auto & kbl1 = jbl1.at(rel_k);
                    const auto & kbl2 = jbl2.at(rel_k);
                    auto & kdbl = jdbl.at(rel_k);

                    for (std::size_t l = k + 1; l < le; ++l) {
                        const auto rel_l = l - k - 1;
                        const auto & el1 = kbl1.at(rel_l);
                        const auto & el2 = kbl2.at(rel_l);
                        const auto qcmp = comp_qt(el1, el2);
                        kdbl.at(rel_l) = qcmp;
                        if (qcmp == Q_COMP::NO_COMP) {
                            continue;
                        }
                        // here we'll just count conflicts as distances,
                        //  so polytomy, compat or same all count as no diff
                        if (qcmp == Q_COMP::CONFLICT_RES) {
                            diff_by_taxon.at(i) += 1;
                            diff_by_taxon.at(j) += 1;
                            diff_by_taxon.at(k) += 1;
                            diff_by_taxon.at(l) += 1;
                            num_diffs += 1;
                        }
                        num_comp += 1;
                        comp_by_taxon.at(i) += 1;
                        comp_by_taxon.at(j) += 1;
                        comp_by_taxon.at(k) += 1;
                        comp_by_taxon.at(l) += 1;                        
                    }
                }
            }
        }
    }
};

void quartet_dist_analysis(const Tree_t & inp_tre1,
                           const Tree_t & inp_tre2) {
    TreeAsUIntSplits tas_1{inp_tre1};
    TreeAsUIntSplits tas_2{inp_tre2};
    if (tas_1.leaf_label_to_ind != tas_2.leaf_label_to_ind) {
        throw OTCError() << "trees must have the same leaf label set.\n";
    }
    AllQuartets t_1_q{tas_1};
    AllQuartets t_2_q{tas_2};
    //t_1_q.write(std::cout);
    QuartDist qdist{t_1_q, t_2_q};
    const auto dc = qdist.get_diff_comp();
    std::cout << dc.first << "\t" << dc.second << "\t" << frac_diff_from_pair(dc) << std::endl;
    throw OTCError() << "Early exit.\n";
    
}


int main(int argc, char *argv[]) {
    try {
        variables_map args = parse_cmd_line(argc, argv);
        string first = args["first"].as<string>();
        string second = args["second"].as<string>();
        ParsingRules rules;
        // 1. Load and process summary tree.
        auto fir_trees = get_trees<Tree_t>(first, rules);
        auto sec_trees = get_trees<Tree_t>(second, rules);
        std::cerr << fir_trees.size() << " trees in " << first << std::endl;
        std::cerr << sec_trees.size() << " trees in " << second << std::endl;
        if (fir_trees.size() != sec_trees.size()) {
            std::cerr << "otc-quartet-distances: Error tree files must have the same number of trees" << std::endl;
            exit(1);
        }
        for (std::size_t ind = 0 ; ind < fir_trees.size(); ++ind) {
            const Tree_t & fir_tree = *(fir_trees.at(ind));
            const Tree_t & sec_tree = *(sec_trees.at(ind));
            quartet_dist_analysis(fir_tree, sec_tree);
        }
        /*
        if (fir_trees.size() < 2) {
            std::cerr << "otc-contrast-astral-runs: Error tree files must have at least 2 trees" << std::endl;
            exit(1);
        }
        const Tree_t & fir_sp_tree = *(fir_trees[0]);
        const TreeAsSplits tas_1(fir_sp_tree);
        const Tree_t & sec_sp_tree = *(sec_trees[0]);
        const TreeAsSplits tas_2(sec_sp_tree);
        if (tas_2.leaf_labels != tas_1.leaf_labels) {
            auto d = set_sym_difference_as_set(tas_2.leaf_labels, tas_1.leaf_labels);
            std::cerr << "otc-contrast-astral-runs: Some labels uniq to one tree:" << std::endl;
            std::cerr << " size of diff set = " << d.size() << "\n";
            for (auto item : d) {
                std::cerr << "  " << item << "\n";
            }
            exit(1);
        }
        const auto & full_leaf_set = tas_1.leaf_labels;
        const auto & to_inf_1 = tas_1.inf_taxset_to_nd;
        const auto & to_inf_2 = tas_2.inf_taxset_to_nd;
        auto all_splits = set_union_as_set(keys(to_inf_1), keys(to_inf_2));
        split_to_summary s2ss;
        for (const auto & sp : all_splits) {
            if (contains(to_inf_1, sp)) {
                if (contains(to_inf_2, sp)) {
                    s2ss.emplace(sp, in_full_status::IN_BOTH);
                } else {
                    s2ss.emplace(sp, in_full_status::IN_ONE);
                }
            } else {
                s2ss.emplace(sp, in_full_status::IN_TWO);
            }
        }
        // check_cs(s2ss);
        unsigned inp_tree_index = 1;
        while (inp_tree_index < fir_trees.size()) {
            const auto & inp_tre1 = *(fir_trees[inp_tree_index]);
            const auto & inp_tre2 = *(sec_trees[inp_tree_index]);
            analyze_inp_tree_pair(tas_1, tas_2, inp_tre1, inp_tre2, s2ss, inp_tree_index, false);
            std::cerr << "tree " << inp_tree_index << std::endl;
            //check_cs(s2ss);
            inp_tree_index++;
        }
        check_cs(s2ss);
        report_summaries(s2ss, "splits-key.json", "split-summary.tsv");
        */


    } catch (std::exception& e) {
        std::cerr << "otc-quartet-distances: Error! " << e.what() << std::endl;
        exit(1);
    }
}
