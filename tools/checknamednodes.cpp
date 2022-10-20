#include <algorithm>
#include <unordered_set>
#include <iterator>
#include <sstream>
#include <stack>
#include <tuple>
#include <deque>
#include <utility>
#include <memory>
#include "otc/otcli.h"
#include "otc/tree_operations.h"

using namespace otc;


static std::string incertaeSedisFilename;

struct RTNodeTrav;
using Node_t = RootedTreeNode<RTNodeTrav>;

struct RTNodeTrav {
    // mutable std::uint32_t trav_enter = UINT32_MAX;
    // mutable std::uint32_t trav_exit = UINT32_MAX;
};

struct RTTreeEmpty {
};

//// Impl
using std::string;
using std::get;
using std::unique_ptr;
using std::vector;

using Tree_t = RootedTree<RTNodeTrav, RTTreeEmpty>;

template<typename N>
OttId get_first_anc_with_id(const N *nd) {
    for (auto anc : iter_anc(*nd)) {
        if (anc->has_ott_id()) {
            return anc->get_ott_id();
        }
    }
    throw OTCError() << "No node on path to root had an OTT ID";
}

template<typename N>
OttId get_first_anc_with_rel_id(const N *nd, const OttIdSet & relevant) {
    for (auto anc : iter_anc(*nd)) {
        if (anc->has_ott_id()) {
            auto ott_id = anc->get_ott_id(); 
            if (contains(relevant, ott_id)) {
                return ott_id;
            }
        }
    }
    throw OTCError() << "No node on path to root had a relevant OTT ID";
}


using pair_ids_sets = std::pair<OttIdSet, OttIdSet>;
using id_2_anc_t = std::map<OttId, OttId>;
template<typename T>
pair_ids_sets get_internal_ids(T & tree,
                        std::map<OttId, OttId> & toNamedAnc) {
    pair_ids_sets pis;
    OttIdSet & tip_ids = pis.first;
    OttIdSet & internal_ids = pis.second;

    for (auto nd : iter_pre(tree)) {
        if (nd->is_tip()) {
            if (!nd->has_ott_id()) {
                throw OTCError() << "tip without OTT ID";
            }
            auto ott_id = nd->get_ott_id();
            tip_ids.insert(ott_id);
            auto first_anc_id = get_first_anc_with_id(nd);
            toNamedAnc[ott_id] = first_anc_id;
        } else {
            if (nd->has_ott_id()) {
                auto ott_id = nd->get_ott_id();
                internal_ids.insert(ott_id);
                if (nd->get_parent() != nullptr) {
                    auto first_anc_id = get_first_anc_with_id(nd);
                    toNamedAnc[ott_id] = first_anc_id;
                }
            }
        }
    }
    return pis;
}

template<typename T>
pair_ids_sets taxo_get_internal_ids(T & tree,
                                    std::map<OttId, OttId> & toNamedAnc,
                                    std::map<OttId, const Node_t *> & id2ndMap,
                                    const OttIdSet & relevant) {
    pair_ids_sets pis;
    OttIdSet & tip_ids = pis.first;
    OttIdSet & internal_ids = pis.second;

    for (auto nd : iter_pre(tree)) {
        if (nd->is_tip()) {
            if (!nd->has_ott_id()) {
                throw OTCError() << "tip without OTT ID";
            }
            auto ott_id = nd->get_ott_id();
            id2ndMap[ott_id] = nd;
            tip_ids.insert(ott_id);
            auto first_anc_id = get_first_anc_with_rel_id(nd, relevant);
            toNamedAnc[ott_id] = first_anc_id;
        } else {
            if (!nd->has_ott_id()) {
                throw OTCError() << "taxonomy internal node without OTT ID";
            }
            auto ott_id = nd->get_ott_id();
            id2ndMap[ott_id] = nd;
            if (contains(relevant, ott_id)) {
                internal_ids.insert(ott_id);
                if (nd->get_parent() != nullptr) {
                    auto first_anc_id = get_first_anc_with_rel_id(nd, relevant);
                    toNamedAnc[ott_id] = first_anc_id;
                }
            }
        }
    }
    return pis;
}

int check_ancs_for_id(const OttId ott_id,
                      OttId root_id,
                      const id_2_anc_t & phylo2anc,
                      const id_2_anc_t & taxo2anc,
                      std::map<OttId, const Node_t *> & tax_id2nd,
                      const OttIdSet & incertae_sedis_ids) {
    if (ott_id == root_id) {
        return 0;
    }
    OttId phy_par_id, tax_par_id;
    try {
        phy_par_id = phylo2anc.at(ott_id);
    } catch(...) {
        throw OTCError() << ott_id << " not in phylo2anc map\n";
    }
    try {
        tax_par_id = taxo2anc.at(ott_id);
    } catch(...) {
        throw OTCError() << ott_id << " not in taxo2anc map\n";
    }
    
    if (phy_par_id == tax_par_id) {
        return 0;
    }
    LOG(DEBUG) << ott_id << " is first a descendant of " << phy_par_id << " in phylogeny, but " << tax_par_id << " in taxonomy.\n";
    
    const Node_t * deepest_inc_sed_nd = nullptr;
    const Node_t * curr_nd = tax_id2nd.at(ott_id);
    OttId nd_id = ott_id;
    while (true) {
        if (nd_id == tax_par_id) {
            break;
        }
        if (contains(incertae_sedis_ids, nd_id)) {
            LOG(DEBUG) << "  setting deepest_inc_sed_nd to " << ott_id;
            deepest_inc_sed_nd = curr_nd; 
        } else {
            LOG(DEBUG) << nd_id << " not incertae_sedis";
        }
        curr_nd = curr_nd->get_parent();
        nd_id = curr_nd->get_ott_id();
    } 
    OttId must_be_anc_id;
    if (deepest_inc_sed_nd != nullptr) {
        must_be_anc_id = deepest_inc_sed_nd->get_parent()->get_ott_id();
        LOG(DEBUG) << "must_be_anc_id = "<< must_be_anc_id;
        OttId phy_spike_id = phy_par_id;
        const Node_t * phy_spike_nd = tax_id2nd.at(phy_spike_id);
        while (true) {
            if (phy_spike_id == must_be_anc_id) {
                LOG(DEBUG) << "must_be_anc_id = " << must_be_anc_id << " encountered";
                return 0;
            }
            LOG(DEBUG) << "phy_spike_id = "<< phy_spike_id << " !=  " << must_be_anc_id;
            if (phy_spike_id == tax_par_id) {
                LOG(DEBUG) << "phy_spike_id = "<< phy_spike_id << " ==  " << tax_par_id << " = tax_par_id. breaking...";
                break;
            }
            phy_spike_nd = phy_spike_nd->get_parent();
            phy_spike_id = phy_spike_nd->get_ott_id();
            if (phy_spike_id == root_id) {
                LOG(DEBUG) << "phy_spike_id = "<< phy_spike_id << " ==  " << root_id << " = root_id, breaking...";
                break;
            }
        }
    } else {
        LOG(DEBUG) << "deepest_inc_sed_nd = nullptr";
    }
    std::cout << ott_id << " is first a descendant of " << phy_par_id << " in phylogeny, but " << tax_par_id << " in taxonomy.\n";
    throw OTCError() << "Doh!";
    return 1;
}

int check_named_nodes(const Tree_t & supertree,
                      const Tree_t & taxonomy,
                      const OttIdSet & incertae_sedis_ids) {
    LOG(DEBUG) << "incertae_sedis_ids.size() = " << incertae_sedis_ids.size();
    int errs = 0;
    id_2_anc_t phylo2anc;
    auto [phylo_tips, phylo_internals] = get_internal_ids(supertree, phylo2anc);
    id_2_anc_t taxo2anc;
    std::map<OttId, const Node_t *> tax_id2nd;
    auto [taxo_tips, taxo_internals] = taxo_get_internal_ids(taxonomy, taxo2anc, tax_id2nd, phylo_internals);
    if (phylo_tips != taxo_tips) {
        auto phylo_extras = set_difference_as_set(phylo_tips, taxo_tips);
        for (auto pe : phylo_extras) {
            errs += 1;
            std::cout << pe << "\textra tip in phylogeny.\n";
        }
        auto taxo_extras = set_difference_as_set(taxo_tips, phylo_tips);
        for (auto te : taxo_extras) {
            errs += 1;
            std::cout << te << "\ttip missing in phylogeny.\n";
        }
    }
    if (phylo_internals != taxo_internals) {
        auto phylo_extras = set_difference_as_set(phylo_internals, taxo_internals);
        for (auto pe : phylo_extras) {
            errs += 1;
            std::cout << pe << "\textra internal in phylogeny.\n";
        }
        auto taxo_extras = set_difference_as_set(taxo_internals, phylo_internals);
        for (auto te : taxo_extras) {
            errs += 1;
            std::cout << te << "\tinternal missing in phylogeny.\n";
        }
    }
    auto root_id = supertree.get_root()->get_ott_id();
    for (auto ott_id : phylo_tips) {
        errs += check_ancs_for_id(ott_id, root_id, phylo2anc, taxo2anc, tax_id2nd, incertae_sedis_ids);
    }
    for (auto ott_id : phylo_internals) {
        errs += check_ancs_for_id(ott_id, root_id, phylo2anc, taxo2anc, tax_id2nd, incertae_sedis_ids);
    }
    return errs;
}


bool handleIncertaeSedis(OTCLI&, const string & arg) {
    incertaeSedisFilename = arg;
    return true;
}


int main(int argc, char *argv[]) {
    OTCLI ot_cli("otc-named-unnamed-nodes",
                "takes a labelled phylogeny and a taxonomy tree (the second tree file)\n"
                "  and a filepath to list of incertae sedis taxa (1 ID per line)."
                "Emits an error message for every taxon that does not have the expected\n"
                "  parent, based on the structure of the taxonomy and the nodes with \n"
                "  OTT IDs in the phylogeny.\n",
                "-iincertae_sedis.txt phylo.tre taxo.tre");
    ot_cli.add_flag('i',
              "Optional list of IDs of tree in the taxonomy that are incertae sedis",
              handleIncertaeSedis,
              true);
    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};
    if (argc < 2) {
        throw OTCError("No trees provided!");
    }
    if (tree_processing_main<Tree_t>(ot_cli, argc, argv, get, nullptr, nullptr, 1)) {
        return 1;
    }
    if (trees.size() != 2) {
        std::cerr << "Supplied " << trees.size() << " trees for checking - that should be 2 trees!";
        return 2;
    }
    OttIdSet incertae_sedis_ids;
    if (!incertaeSedisFilename.empty()) {
        std::ifstream incert_sed_id_file(incertaeSedisFilename.c_str());
        if (!incert_sed_id_file.good()) {
            std::cerr << "Could not open the incertae sedis ID file: " << incertaeSedisFilename << "\n";
            return 1;
        }
        string line;
        while (std::getline(incert_sed_id_file, line)) {
            char* temp;
            long raw_ott_id = std::strtoul(line.c_str(), &temp, 10);
            if (*temp != '\0' && *temp != '\n') {
                std::cerr << "Expecting just numbers and newlines in incertae sedis file found: " << line << "\n";
                return 1;
            }
            incertae_sedis_ids.insert(check_ott_id_size(raw_ott_id));
        }
    }
    auto & supertree = *(trees.at(0));
    auto & taxonomy = *(trees.at(1));
    auto rc = 1;
    try {
        rc = check_named_nodes(supertree, taxonomy, incertae_sedis_ids);
    } catch (OTCError & err) {
        std::cerr << "Exception: " << err.what() << std::endl;
        rc = 1;
    }
    return rc;
}
