#ifndef OTCETERA_ALL_QUARTETS_H
#define OTCETERA_ALL_QUARTETS_H
#include <ostream>
#include <string>
#include <vector>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/tree_as_split_set.h"

namespace otc {

template <typename T> class QuartDist;
enum QUARTET_TYPE {UNKNOWN = 0,
                   POLYTOMY = 1,
                   ONE_TWO = 2,
                   ONE_THREE = 3,
                   ONE_FOUR = 4,
                   NOT_Q = 5
                   };

void write_qt(std::ostream & out, const QUARTET_TYPE & qt);

template<typename T>
std::vector<T> gen_by_fourth(std::size_t num_tax,
							 std::size_t third_index,
							 const T & def);
template<typename T>
std::vector<std::vector<T> > gen_by_third(std::size_t num_tax,
										  std::size_t sec_index,
										  const T & def );

template<typename T>
std::vector<std::vector<std::vector<T> > > gen_by_sec(std::size_t num_tax,
													  std::size_t first_ind,
													  const T & def);


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
gen_by_third(std::size_t num_tax, std::size_t sec_index, const T & def ) {
    std::vector<std::vector<T> > ret;
    const std::size_t min_real_ind = sec_index + 1;
    assert(min_real_ind < num_tax);
    // N - 1 (for 2 more indices) - 1 (for the zero-based ind)
    const std::size_t max_real_ind = num_tax - 2;
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
    std::vector<std::vector<std::vector<T> > > ret;
    const std::size_t min_real_ind = first_ind + 1;
    assert(min_real_ind < num_tax);
    // N - 2 (for 2 more indices) - 1 (for the zero-based ind)
    const std::size_t max_real_ind = num_tax - 3;
    assert(max_real_ind >= min_real_ind);
    const std::size_t my_size = 1 + max_real_ind - min_real_ind;
    ret.reserve(my_size);
    for (std::size_t my_real_ind = min_real_ind; my_real_ind <=  max_real_ind; ++my_real_ind) {
        ret.push_back(gen_by_third(num_tax, my_real_ind, def));
    }
    return ret;
}


template <typename T>
class AllQuartets {
    using off_by_3_v = std::vector<QUARTET_TYPE> ;
    using off_by_2_v = std::vector<off_by_3_v>;
    using off_by_1_v = std::vector<off_by_2_v>;
    using top_level = std::vector<off_by_1_v>;

    top_level by_lowest;

    public:
    using uint_set = std::set<std::size_t>;
    using node_t = typename T::node_type;

    std::size_t num_tips;

    AllQuartets(const GenTreeAsUIntSplits<T> & tas) {
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
            register_polytomy(par_nd, outgroup, tas);
        }
    }

    void register_polytomy(const node_t * par_nd, 
                           const uint_set & outgroup,
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

    friend class QuartDist<T>;
};


} // namespace otc

#endif
