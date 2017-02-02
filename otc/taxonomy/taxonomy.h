#ifndef OTC_TAXONOMY_TAXONOMY_H
#define OTC_TAXONOMY_TAXONOMY_H
// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <boost/utility/string_ref.hpp>
#include <bitset>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"

// 1. Write a parser to read the lines faster
// 2. Avoid memory allocation -- by mmapping the taxonomy file?
// 3. Convert the flags into a bitmask
// 4. Should the Rank be a converted to an integer?
// 5. Can we assign OTT IDs to internal nodes of a tree while accounting for Incertae Sedis taxa?
// * What are the triplet-inference rules for the Incertae Sedis problem?

namespace otc {

struct taxonomy_record {
    std::string line;
    long id = 0;
    long parent_id = 0;
    int parent_index = 0;
    boost::string_ref name;
    boost::string_ref rank;
    boost::string_ref sourceinfo;
    boost::string_ref uniqname;
    std::bitset<32> flags;
    int depth = 0;
    int out_degree = 0;
    taxonomy_record(taxonomy_record&& tr) = default;
    explicit taxonomy_record(const std::string& line);
    std::vector<std::string> sourceinfoAsVec() const {
        std::vector<std::string> mt;
        if (sourceinfo.empty()) {
            return mt;
        }
        std::string sis = std::string(sourceinfo);
        auto aslist = split_string(sis, ',');
        mt.assign(aslist.begin(), aslist.end());
        return mt;
    }
};

struct Taxonomy: public std::vector<taxonomy_record> {
    std::unordered_map<long,int> index;
    std::unordered_map<long,long> forwards;
    long keep_root;
    std::bitset<32> cleaning_flags;
    std::string path;
    std::string version;
    template <typename Tree_t> std::unique_ptr<Tree_t> getTree(std::function<std::string(const taxonomy_record&)>) const;

          taxonomy_record& record_from_id(long id);
    const taxonomy_record& record_from_id(long id) const;
    
    /// Index of the root record
    constexpr static int root_index() {
        return 0;
    }

    /// Given an OTT ID, forward if necessary and return the current ID.  Returns -1 if not found.
    long map(long id) const;

    /// Write out a taxonomy to directory dirname
    void write(const std::string& dirname);

    /// Load the taxonomy from directory dir, and apply cleaning flags cf, and keep subtree below kr
    Taxonomy(const std::string& dir, std::bitset<32> cf = std::bitset<32>(), long kr = -1);
};

class RTTaxNodeData {
    public:
    const taxonomy_record * taxonomy_line = nullptr;
};

template <typename Node_t>
void populateNodeFromTaxonomyRecord(Node_t & nd,
                                    const taxonomy_record & line,
                                    std::function<std::string(const taxonomy_record&)> getName);


// default behavior is to set ID and Name from line
template <typename Node_t>
inline void populateNodeFromTaxonomyRecord(Node_t & nd,
                                           const taxonomy_record & line,
                                           std::function<std::string(const taxonomy_record&)> getName) {
    nd.setOttId(line.id);
    nd.setName(getName(line));    
}

// default behavior is to set ID and Name from line
template <>
inline void populateNodeFromTaxonomyRecord(RootedTreeNode<RTTaxNodeData> & nd,
                                           const taxonomy_record & line,
                                           std::function<std::string(const taxonomy_record&)>) {
    nd.setOttId(line.id);
    nd.getData().taxonomy_line = &line;    
}


template <typename Tree_t>
std::unique_ptr<Tree_t> Taxonomy::getTree(std::function<std::string(const taxonomy_record&)> getName) const {
    const auto& taxonomy = *this;
    std::unique_ptr<Tree_t> tree(new Tree_t);
    vector<typename Tree_t::node_type*> node_ptr(size(), nullptr);
    for(auto i = 0U; i < taxonomy.size() ; i++) {
        const auto& line = taxonomy[i];
        // Make the tree
        typename Tree_t::node_type* nd = nullptr;
        if (i==0) {
            nd = tree->createRoot();
        } else {
            auto parent_nd = node_ptr[line.parent_index];
            nd = tree->createChild(parent_nd);
        }
        populateNodeFromTaxonomyRecord(*nd, line, getName);
        node_ptr[i] = nd;
    }
    return tree;
}

long root_ott_id_from_file(const std::string& filename);
std::string get_taxonomy_dir(const boost::program_options::variables_map& args);
Taxonomy load_taxonomy(const boost::program_options::variables_map& args);

} // namespace
#endif
