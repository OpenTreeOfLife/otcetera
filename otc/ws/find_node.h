#ifndef FIND_NODE_H
#define FIND_NODE_H

#include "tolws.h"
#include <tuple>

namespace otc {

// BEGIN ******* data TaxonToSynth = TaxonPruned | TaxonBroken | TaxonFound ********** //

struct TaxonPruned { };

struct TaxonBroken
{
    const SumTreeNode_t* mrca = nullptr;
    const SumTreeNode_t* node() const {return mrca;}
};

struct TaxonFound
{
    const SumTreeNode_t* found = nullptr;
    const SumTreeNode_t* node() const {return found;}
};

struct TaxonToSynth: public std::variant<TaxonPruned, TaxonBroken, TaxonFound>
{
    const SumTreeNode_t* node() const
    {
        switch (index())
        {
        case 1: return std::get<TaxonBroken>(*this).node(); break;
        case 2: return std::get<TaxonFound>(*this).node(); break;
        default: return nullptr; break;
        }
        std::abort();
    }
    bool ok() const {return node();}

    bool pruned() const {return index() == 0;}
    bool broken() const {return index() == 1;}
    bool found() const {return index() == 2;}

    using std::variant<TaxonPruned, TaxonBroken, TaxonFound>::variant;
};

// BEGIN       data OTTNameToSynth = BadID | InvalidID OttID | ValidID ID (Maybe ID) TaxonToSynth      //

// NOTE: For too-big, we should probably directly check if we can read the result into a long.
// See https://stackoverflow.com/questions/194465/how-to-parse-a-string-to-an-int-in-c
struct BadID { };

struct InvalidID
{
    OttId id;
};

struct ValidID
{
    OttId id;
    std::optional<OttId> forwarded_from;
    TaxonToSynth to_synth;
    const SumTreeNode_t* node() const
    {
        return to_synth.node();
    }
    bool ok() const { return node();}
    bool pruned() const {return to_synth.pruned();}
    bool broken() const {return to_synth.broken();}
    bool found() const {return to_synth.found();}
};

// BEGIN NameToSynth = NoMatchName | OTTNameToSynth | MRCANameToSynth

struct NoMatchName { };

struct OTTNameToSynth : public std::variant<BadID,InvalidID,ValidID>
{
    const SumTreeNode_t* node() const
    {
        if (index() == 2)
            return std::get<ValidID>(*this).node();
        else
            return nullptr;
    }
    bool ok() const {return node();}

    bool invalid() const {return index() == 1;}
    bool pruned() const {return index() == 2 and std::get<ValidID>(*this).pruned();}
    bool broken() const {return index() == 2 and std::get<ValidID>(*this).broken();}
    bool found() const {return index() == 2 and std::get<ValidID>(*this).found();}

    using std::variant<BadID,InvalidID,ValidID>::variant;
};


struct MRCANameToSynth
{
    OTTNameToSynth node1;
    OTTNameToSynth node2;
    const SumTreeNode_t* mrca = nullptr;
    const SumTreeNode_t* node() const {return mrca;}
    bool ok() const {return node();}
};

struct NameToSynth: public std::variant<NoMatchName, OTTNameToSynth, MRCANameToSynth>
{
    const SumTreeNode_t* node() const
    {
        switch (index())
        {
        case 1: return std::get<OTTNameToSynth>(*this).node(); break;
        case 2: return std::get<MRCANameToSynth>(*this).node(); break;
        default: return nullptr; break;
        }
        std::abort();
    }
    bool ok() const { return node();}

    const OTTNameToSynth& ott_name_lookup() const {return std::get<OTTNameToSynth>(*this);}

    bool is_ott_name() const {return index() == 1;}
    bool is_mrca_name() const {return index() == 2;}
    bool invalid() const {return is_ott_name() and ott_name_lookup().invalid();}
    bool pruned() const {return is_ott_name() and ott_name_lookup().pruned();}
    bool broken() const {return is_ott_name() and ott_name_lookup().broken();}
    bool found() const {return is_ott_name() and ott_name_lookup().found();}

    using std::variant<NoMatchName,OTTNameToSynth,MRCANameToSynth>::variant;
};

std::optional<OttId> is_ott_id(const std::string& node_id);
TaxonToSynth find_node_by_valid_ottid(const SummaryTree_t & tree, OttId id, const std::string& node_id);
OTTNameToSynth find_node_by_ottid_str(const SummaryTree_t & tree, const RichTaxonomy& taxonomy, const std::string & node_id);
MRCANameToSynth find_node_by_mrca_str(const SummaryTree_t & tree, const RichTaxonomy& taxonomy, const std::string & node_id);
NameToSynth find_node_by_id_str(const SummaryTree_t & tree, const RichTaxonomy&, const std::string & node_id);
NameToSynth find_required_node_by_id_str(const SummaryTree_t & tree, const RichTaxonomy& taxonomy, const std::string & node_id);

std::tuple<std::vector<const SumTreeNode_t*>,nlohmann::json>
find_nodes_for_id_strings(const RichTaxonomy& taxonomy, const SummaryTree_t* tree_ptr, const std::vector<std::string>& node_ids,
                          bool fail_broken = false);
std::string find_node_failure_reason(const NameToSynth& result, bool fail_broken = false);
}
#endif
