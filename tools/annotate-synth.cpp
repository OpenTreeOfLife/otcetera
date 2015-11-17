#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include <tuple>
#include <sstream>
#include <cstring>
#include <unordered_map>

using namespace otc;

using std::string;

struct RTNodeDepth
{
    int depth = 0;
};


using Tree_t = RootedTree<RTNodeDepth, RTreeNoData>;
using node_t = Tree_t::node_type;

bool showJSON = false;

/// Stat Calc declarations
enum NDSB {
    ROOT_BIT = 0x01,
    OUTDEGREE_ONE_BIT = 0x02,
    LEAF_BIT = 0x04,
    DISPLAYED_BIT = 0x08,
    COULD_RESOLVE_BIT = 0x10,
    INCOMPATIBLE_BIT = 0x20,
    ANC_ALL_BIT = 0x40,
    END_BIT = 0x80
};
enum NDSE {
    ROOT_NODE = NDSB::ROOT_BIT | NDSB::ANC_ALL_BIT , // root of graph and has out-degree>1
    FIRST_FORK = NDSB::ANC_ALL_BIT , // root of phylo info, but has redundant parent
    LEAF_NODE = NDSB::LEAF_BIT,
    FORKING_DISPLAYED = NDSB::DISPLAYED_BIT,
    FORKING_COULD_RESOLVE = NDSB::COULD_RESOLVE_BIT,
    FORKING_INCOMPATIBLE = NDSB::INCOMPATIBLE_BIT,
    REDUNDANT_DISPLAYED = NDSB::DISPLAYED_BIT | NDSB::OUTDEGREE_ONE_BIT,
    REDUNDANT_COULD_RESOLVE = NDSB::COULD_RESOLVE_BIT | NDSB::OUTDEGREE_ONE_BIT,
    REDUNDANT_INCOMPATIBLE = NDSB::INCOMPATIBLE_BIT | NDSB::OUTDEGREE_ONE_BIT,
    REDUNDANT_TERMINAL = NDSB::LEAF_BIT | NDSB::OUTDEGREE_ONE_BIT,
    // is and anc of FIRST_FORK, but not root of graph
    REDUNDANT_ROOT_ANC = NDSB::ANC_ALL_BIT | NDSB::OUTDEGREE_ONE_BIT,
    // Is root of graph, but anc of FIRST_FORK
    ROOT_REDUNDANT_ROOT_ANC = NDSB::ROOT_BIT | NDSB::ANC_ALL_BIT | NDSB::OUTDEGREE_ONE_BIT,
    //ROOT is a TIP
    DOT_TREE = NDSB::ROOT_BIT | NDSB::ANC_ALL_BIT | NDSB::LEAF_BIT,
    REDUNDANT_LINE_TREE = NDSB::ANC_ALL_BIT | NDSB::LEAF_BIT | NDSB::OUTDEGREE_ONE_BIT,
    ROOT_REDUNDANT_LINE_TREE = NDSB::ANC_ALL_BIT | NDSB::LEAF_BIT | NDSB::OUTDEGREE_ONE_BIT,
    END_VALUE = NDSB::END_BIT
};

inline NDSE operator|(NDSE f, NDSE s) {
    return static_cast<NDSE>(static_cast<int>(f) | static_cast<int>(s));
}
inline NDSE operator|(NDSB f, NDSE s) {
    return static_cast<NDSE>(static_cast<int>(f) | static_cast<int>(s));
}
inline NDSE operator|(NDSE f, NDSB s) {
    return static_cast<NDSE>(static_cast<int>(f) | static_cast<int>(s));
}

long n_leaves(const node_t* nd)
{
    long count = 0;
    for(auto nd2: iter_post_n_const(*nd))
        if (nd2->isTip())
            count++;
    return count;
}

string study_from_tree_name(const string& name)
{
    const char* start = strrchr(name.c_str(),' ');
    start++;
    assert(start);

    const char* end = strrchr(name.c_str(),'_');
    return name.substr(start - name.c_str(), end-start);
}


string tree_in_study_from_tree_name(const string& name)
{
    const char* start = strrchr(name.c_str(),'_');
    start++;
    assert(start);

    const char* end = strrchr(name.c_str(),'.');
    return name.substr(start - name.c_str(), end-start);
}

string getNodeName(const string& name)
{
    const char* start = name.c_str();
    const char* end = name.c_str() + name.size();
    while(strchr(" \t_",*start) and start < end)
        start++;
    while(strchr(" \t_",*(end-1)) and start < end)
        end--;
    if (start == end)
        throw OTCError()<<"Node name '"<<name<<"' contracted to nothing!";
    return name.substr(start - name.c_str(), end-start);
}

string quote(const string& s)
{
    return '"'+s+'"';
};

std::map<NDSE, std::size_t> doStatCalc(const Tree_t & summaryTree,
                                       const Tree_t & inpTree,
                                       std::map<const Tree_t::node_type *, NDSE> * node2Classification,
                                       std::unordered_multimap<string,string> * support,
                                       std::unordered_multimap<string,string> * conflict,
                                       bool isTaxoComp) {
}
/// end Stat Calc impl
/// Stat report decl
void writeHeader(std::ostream &out);
void writeRow(std::ostream &out, std::map<NDSE, std::size_t> & m, const std::string & label);
std::string explainOutput();
/// End Stat report decl
/// Stat report impl
std::string explainOutput() {
    std::ostringstream o;
    o << "Writes tab-separated output.\n";
    o << "Each row reports the number of internal nodes of the input tree that fall into each category.\n";
    o << "The final row shows the totals.\n";
    o << "The two \"axes\" that the statistics explore are support and out-degree.\n";
    o << "Columns starting with \"F\" are \"forking\" internal nodes with out-degree > 1.\n";
    o << "Columns starting with \"R\" are \"redundant\" internal nodes with out-degree = 1.\n";
    o << "A \"D\" suffix to a column header means that the node is displayed by the summary tree.\n";
    o << "A \"CR\" suffix means that the node is could resolve a polytomy in the summary tree (so the\n";
    o << "    summary tree is not unambiguously in conflict in the node).\n";
    o << "An \"I\" suffix to a column header means that the node is incompatible with every resolution ofthe summary tree.\n";
    o << "For the redundant nodes, the report indicates the conflict status of their closest non-redundant descendant.\n";
    o << "A redundant node can also be marked \"T\" (for \"trivial\")if it is an ancestor of only 1 leaf or of the root.\n";
    o << "The \"F\" and \"R\" column are just the sums for forking and redundant entries.\n";
    o << "The \"label\" shows the tree name or \"Total of # trees\" for the global sum\n";
    o << "The ordering of the rows is the input order. For columns it is:\n";
    o <<  "FD FCR FI F RD RCR RI RT R label";
    return o.str();
}

void writeHeader(std::ostream &out) {
    out << "FD" << '\t'
        << "FCR" << '\t'
        << "FI" << '\t'
        << "F" << '\t'
        << "RD" << '\t'
        << "RCR"  << '\t'
        << "RI" << '\t'
        << "RT" << '\t'
        << "R" << '\t'
        << "label" << '\n';
}


void writeRow(std::ostream &out,
              std::map<NDSE, std::size_t> & m,
              const std::string & label) {
    const auto f = m[NDSE::FORKING_DISPLAYED] + m[NDSE::FORKING_COULD_RESOLVE] + m[NDSE::FORKING_INCOMPATIBLE];
    const auto rt = m[NDSE::REDUNDANT_TERMINAL] + m[NDSE::REDUNDANT_ROOT_ANC];
    const auto r = m[NDSE::REDUNDANT_DISPLAYED] + m[NDSE::REDUNDANT_COULD_RESOLVE] + m[NDSE::REDUNDANT_INCOMPATIBLE] + rt;
    out << m[NDSE::FORKING_DISPLAYED] << '\t'
        << m[NDSE::FORKING_COULD_RESOLVE]  << '\t'
        << m[NDSE::FORKING_INCOMPATIBLE] << '\t'
        << f << '\t'
        << m[NDSE::REDUNDANT_DISPLAYED] << '\t'
        << m[NDSE::REDUNDANT_COULD_RESOLVE] << '\t'
        << m[NDSE::REDUNDANT_INCOMPATIBLE] << '\t'
        << rt << '\t'
        << r  << '\t'
        << label << '\n';
}

struct DisplayedStatsState : public TaxonomyDependentTreeProcessor<Tree_t> {
    std::unique_ptr<Tree_t> summaryTree;
    std::map<NDSE, std::size_t> totals;
    std::unordered_multimap<string,string> support;
    std::unordered_multimap<string,string> conflict;
    int numErrors = 0;
    bool treatTaxonomyAsLastTree = false;
    bool headerEmitted = false;
    int numTrees = 0;
    virtual ~DisplayedStatsState(){}

    bool summarize(OTCLI &otCLI) override {
        if (treatTaxonomyAsLastTree) {
            statsForNextTree(otCLI, *taxonomy, true);
        }
        if (not showJSON)
        {
            const std::string label = std::string("Total of ") + std::to_string(numTrees) + std::string(" trees");
            writeNextRow(otCLI.out, totals, label);
        }
        else
        {
            std::cout<<"  \"nodes\": {\n";
            for(auto nd: iter_post_const(*summaryTree))
            {
                string name = nd->getName();
                if (nd->hasOttId())
                    name = "ott" + std::to_string(nd->getOttId());

                int sc = support.count(name);
                int cc = conflict.count(name);
                if (sc + cc == 0) continue;

                std::cout<<"    "<<quote(name)<<": { \n";
                if (sc)
                {
                    std::cout<<"      \"supported-by\": [ ";
                    auto range = support.equal_range(name);
                    auto start = range.first;
                    auto end   = range.second;
                    for (auto iter = start; iter!=end; ++iter)
                    {
                        if (iter != start) std::cout<<"                        ";
                        std::cout<<iter->second;
                        auto next = iter; ++next;
                        if (next != end)
                            std::cout<<",\n";
                    }
                    std::cout<<" ]";
                    if (cc > 0)
                        std::cout<<",";
                    std::cout<<"\n";
                }
                if (cc)
                {
                    std::cout<<"      \"conflicts-with\": [ ";
                    auto range = conflict.equal_range(name);
                    auto start = range.first;
                    auto end   = range.second;
                    for (auto iter = start; iter!=end; ++iter)
                    {
                        if (iter != start) std::cout<<"                        ";
                        std::cout<<iter->second;
                        auto next = iter; ++next;
                        if (next != end)
                            std::cout<<",\n";
                    }
                    std::cout<<" ]\n";
                }
                std::cout<<"    }\n";
            }
            std::cout<<"  }\n";
        }
        return true;
    }

    void writeNextRow(std::ostream &out,
                      std::map<NDSE, std::size_t> & m,
                      const std::string & label) {
        if (!headerEmitted) {
            writeHeader(out);
            headerEmitted = true;
        }
        writeRow(out, m, label);
    }

    void statsForNextTree(OTCLI & otCLI, const Tree_t & tree, bool isTaxoComp) {
//        auto c = doStatCalc(*summaryTree, tree, nullptr, showJSON?(&support):nullptr, showJSON?(&conflict):nullptr, isTaxoComp);
//        if (not showJSON) writeNextRow(otCLI.out, c, tree.getName());
//        for (const auto & p : c) {
//            totals[p.first] += p.second;
//        }
//        numTrees += 1;
    }

    virtual bool processTaxonomyTree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<Tree_t>::processTaxonomyTree(otCLI);
        otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
        otCLI.getParsingRules().requireOttIds = false;
        // now we get a little cute and reprocess the taxonomy desIds so that they 
        // exclude internals. So that when we expand source trees, we expand just
        // to the taxonomy's leaf set
        return true;
    }

    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) override {
        assert(taxonomy != nullptr);
        if (summaryTree == nullptr) {
            summaryTree = std::move(tree);
            return true;
        }
        requireTipsToBeMappedToTerminalTaxa(*tree, *taxonomy);
        statsForNextTree(otCLI, *tree, false);
        return true;
    }

};


bool handleCountTaxonomy(OTCLI & otCLI, const std::string &) {
    DisplayedStatsState * proc = static_cast<DisplayedStatsState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->treatTaxonomyAsLastTree = true;
    return true;
}

bool handleJSON(OTCLI & otCLI, const std::string &) {
    showJSON = true;
    return true;
}

int main(int argc, char *argv[]) {
    std::string explanation{"takes at least 2 newick file paths: a taxonomy,  a full supertree, and some number of input trees.\n"};
    explanation += explainOutput();
    OTCLI otCLI("otc-displayed-stats",
                explanation.c_str(),
                "synth.tre inp1.tre inp2.tre ...");
    DisplayedStatsState proc;
    otCLI.addFlag('x',
                  "Automatically treat the taxonomy as an input in terms of supporting groups",
                  handleCountTaxonomy,
                  false);
    otCLI.addFlag('j',
                  "Output JSON for node support, instead of displaying statistics.",
                  handleJSON,
                  false);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
