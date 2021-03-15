#ifndef OTCETERA_ALL_TRIPLETS_H
#define OTCETERA_ALL_TRIPLETS_H
#include <ostream>
#include <string>
#include <vector>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/tree_as_split_set.h"

namespace otc {

template <typename T> class TripletDist;
enum TRIPLET_TYPE {UNKNOWN = 0,
                   POLYTOMY = 1,
                   ONE_TWO = 2,
                   ONE_THREE = 3,
                   TWO_THREE = 4,
                   NOT_Q = 5
                   };

void write_tt(std::ostream & out, const TRIPLET_TYPE & qt);

template<typename T>
std::vector<T> gen_tm_by_third(std::size_t num_tax,
							   std::size_t third_index,
							   const T & def);
template<typename T>
std::vector<std::vector<T> > gen_tm_by_sec(std::size_t num_tax,
										  std::size_t sec_index,
										  const T & def );



inline void write_tt(std::ostream & out, const TRIPLET_TYPE & qt) {
    switch(qt) {
        case TRIPLET_TYPE::UNKNOWN : out << "?" ; break;
        case TRIPLET_TYPE::POLYTOMY : out << "*" ; break;
        case TRIPLET_TYPE::NOT_Q : out << "X" ; break;
        case TRIPLET_TYPE::ONE_TWO : out << "12" ; break;
        case TRIPLET_TYPE::ONE_THREE : out << "13" ; break;
        case TRIPLET_TYPE::TWO_THREE : out << "23" ; break;
    }
}

template<typename T>
inline std::vector<T> 
gen_tm_by_third(std::size_t num_tax, std::size_t third_index, const T & def) {
    const std::size_t min_real_ind = third_index + 1;
    assert(min_real_ind < num_tax);
    // N - 1 (for the zero-based ind)
    const std::size_t max_real_ind = num_tax - 1;
    assert(max_real_ind >= min_real_ind);
    const std::size_t my_size = 1 + max_real_ind - min_real_ind;
    std::vector<T> ret{my_size, def};
    return ret;
}

template<typename T>
inline std::vector<std::vector<T> >
gen_tm_by_sec(std::size_t num_tax, std::size_t sec_index, const T & def ) {
    std::vector<std::vector<T> > ret;
    const std::size_t min_real_ind = sec_index + 1;
    assert(min_real_ind < num_tax);
    // N - 1 (for 2 more indices) - 1 (for the zero-based ind)
    const std::size_t max_real_ind = num_tax - 2;
    assert(max_real_ind >= min_real_ind);
    const std::size_t my_size = 1 + max_real_ind - min_real_ind;
    ret.reserve(my_size);
    for (std::size_t my_real_ind = min_real_ind; my_real_ind <=  max_real_ind; ++my_real_ind) {
        ret.push_back(gen_tm_by_third(num_tax, my_real_ind, def));
    }
    return ret;
}

template <typename T>
class AllTriplets {
    using off_by_2_v = std::vector<TRIPLET_TYPE> ;
    using off_by_1_v = std::vector<off_by_2_v>;
    using top_level = std::vector<off_by_1_v>;

    top_level by_lowest;

    public:
    using uint_set = std::set<std::size_t>;
    using node_t = typename T::node_type;

    std::size_t num_tips;

    AllTriplets(const GenTreeAsUIntSplits<T> & tas) {
        num_tips = tas.ind_to_nd.size();
        if (num_tips < 3) {
            return;
        }
        std::size_t n_r = static_cast<std::size_t>(static_cast<int>(num_tips) - 2);
        by_lowest.reserve(n_r);
        const TRIPLET_TYPE def = TRIPLET_TYPE::UNKNOWN;
        for (std::size_t row_n = 0; row_n < n_r; ++row_n) {
            by_lowest.push_back(gen_tm_by_sec<TRIPLET_TYPE>(num_tips, row_n, def));
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
            write_tm_by_sec(out, row, i++);
        }
    }

    private:

    void register_nd(const node_t * par_nd, 
                     const uint_set & outgroup,
                     const GenTreeAsUIntSplits<T> & tas) {
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
            register_polytomy(par_nd, tas);
        }
    }

    void register_polytomy(const node_t * par_nd, 
                           const GenTreeAsUIntSplits<T> & tas) {
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
                    this->register_poly_out(cc_ind_set, cf_ind_set, cs_ind_set);
                    const node_t * curr_t_sib = curr_s_sib->get_next_sib();
                    curr_s_sib = curr_s_sib->get_next_sib();
                }
                curr_f_sib = curr_f_sib->get_next_sib();
            }
            curr_ch = curr_ch->get_next_sib();
        }
    }

    void register_poly_out(const uint_set & f_ind_set, 
                           const uint_set & s_ind_set,
                       	   const uint_set & t_ind_set) {
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
                    this->register_sorted(TRIPLET_TYPE::POLYTOMY, fst_small, fst_mid, fst_large);
                }
            }
        }   
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
                    register_triplet(in_small, in_large, *oi_it);
                }
            }
        }
    }

    void register_triplet(std::size_t in_small, std::size_t in_large, std::size_t out_small) {
        assert(in_small < in_large);
        if (in_small < out_small) {
            if (in_large < out_small) {
                this->register_sorted(TRIPLET_TYPE::ONE_TWO, in_small, in_large, out_small);
            } else {
                this->register_sorted(TRIPLET_TYPE::ONE_THREE, in_small, out_small, in_large);
            }
        } else {
            assert(out_small < in_small);
            this->register_sorted(TRIPLET_TYPE::TWO_THREE, out_small, in_small, in_large);
        }
    }

    void register_sorted(TRIPLET_TYPE tt, std::size_t fir_ind, std::size_t sec_ind, std::size_t thi_ind) {
        assert(thi_ind > sec_ind);    
        std::size_t rel_thi = thi_ind - sec_ind - 1;
        assert(sec_ind > fir_ind);
        std::size_t rel_sec = sec_ind - fir_ind - 1;
        try {
            this->by_lowest.at(fir_ind).at(rel_sec).at(rel_thi) = tt;
        } catch (...) {
            std::cerr << "Error registering ";
            write_tt(std::cerr, tt);
            std::cerr << " at (" << fir_ind << ", " << sec_ind << ", " << thi_ind << ")" << std::endl;
            throw;
        }
    }
    

    void write_tm_by_sec(std::ostream & out, const off_by_1_v & by_sec, std::size_t first) const {
        std::size_t i = first + 1;
        // out << by_sec.size() << " 2nd-level rows\n";
        for (const auto & row : by_sec) {
            write_tm_by_third(out, row, first, i++);
        }
    }

    void write_tm_by_third(std::ostream & out, const off_by_2_v & by_third, std::size_t first, std::size_t sec) const {
        std::size_t i = sec + 1;
        for (const TRIPLET_TYPE el : by_third) {
            if (i > sec + 1) {
                out << " ";
            }
            out << "rt(" << first << ", " << sec << ", " << i << ")=";
            write_tt(out, el);
            i++;
        }
        out << '\n';
    }

    friend class TripletDist<T>;
};


} // namespace otc

#endif
