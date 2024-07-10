#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <bitset>
#include <regex>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "otc/config_file.h"

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::unique_ptr;
using std::set;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;

variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;
    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("taxonomy", value<string>(),"Filename for the taxonomy")
        ;

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("config,c",value<string>(),"Config file containing flags to filter")
        ("clean",value<string>(),"Comma-separated string of flags to filter")
        ("root,r", value<OttId>(), "OTT id of root node of subtree to keep")
        ("xroot,x", value<OttId>(), "OTT id of root node of subtree to keep and detach from parent by writing empty parent id (this is only relevant when the --write-taxonomy option in effect)")
        ;

    options_description selection("Selection options");
    selection.add_options()
    ("cull-flags",value<string>(),"Show records with none of these flags")
    ("any-flags",value<string>(),"Show records with one of these flags")
    ("all-flags",value<string>(),"Show records with all of these flags")
    ("in-tree",value<string>(),"Show records from OTT ids in tree")
    ("in-file",value<string>(),"Show records of integer ids in file");

    options_description output("Output options");
    output.add_options()
        ("show-root,R","Show the ottid of the root node")
        ("find,S",value<string>(),"Show taxa whose names match regex <arg>")
        ("id,I",value<string>(),"Show taxon with OTT ID <arg>")
        ("degree,D",value<long>(),"Show out the degree of node <arg>")
        ("extinct-to-incert,E","Adds an incertae_sedis flag to every extinct taxa (use with write-taxonomy)")
        ("children,C",value<long>(),"Show the children of node <arg>")
        ("parent,P",value<OttId>(),"Show the parent taxon of node <arg>")
        ("report-dist-to-root","Report the number of nodes between from the root to each OTT ID (the root reports 1)")
        ("high-degree-nodes",value<int>(),"Show the top <arg> high-degree nodes")
        ("write-tree,T","Write out the result as a tree")
        ("write-taxonomy",value<string>(),"Write out the result as a taxonomy to directory 'arg'")
        ("name,N", value<OttId>(), "Return name of the given ID")
        ("uniqname,U", value<OttId>(), "Return unique name for the given ID")
        ("report-lost-taxa",value<string>(), "A tree to report missing taxa for")
        ("version,V","Taxonomy version")
        ;

    options_description formatting("Formatting options");
    formatting.add_options()
    ("format",value<string>()->default_value("ott%I\t'%U'\tflags=%F"),"Form of line to write for each taxonomy record");

    options_description visible;
    visible.add(taxonomy).add(selection).add(output).add(formatting).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("taxonomy", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-taxonomy-parser <taxonomy-dir> [OPTIONS]\n"
                                                    "Read a taxonomy, clean it, and then make a tree or some other operation.",
                                                    visible, invisible, p);

    return vm;
}

void report_lost_taxa(const Taxonomy& taxonomy, const string& filename) {
    vector<unique_ptr<Tree_t>> trees;
    std::function<bool(unique_ptr<Tree_t>)> a = [&](unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;};
    ParsingRules rules;
    rules.require_ott_ids = false;
    otc::process_trees(filename,rules,a);//[&](unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;});
    const auto& T =  trees[0];
    std::unordered_map<long, const Tree_t::node_type*> ottid_to_node;
    for(auto nd: iter_pre_const(*T)) {
        if (nd->has_ott_id()) {
            ottid_to_node[nd->get_ott_id()] = nd;
        }
    }
    vector<const TaxonomyRecord*> records;
    for(const auto& rec: taxonomy) {
        records.push_back(&rec);
    }
    std::sort(records.begin(), records.end(), [](const auto& a, const auto& b) {return a->depth < b->depth;});
    for(const auto& rec: records){
        if (not ottid_to_node.count(rec->id)){
            std::cout << "depth=" << rec->depth << "   id=" << rec->id << "   uniqname='" << rec->uniqname << "'\n";
        }
    }
}

void show_rec(const TaxonomyRecord& rec) {
    std::cout << rec.id << "   '" << rec.uniqname << "'   '" << rec.rank << "'   depth = " << rec.depth << "   out-degree = " << rec.out_degree << "    flags = " << flags_to_string(rec.flags) << "\n";
}

vector<OttId> get_ids_from_stream(std::istream& file) {
    vector<OttId> ids;
    long i;
    while (file >> i) {
        ids.push_back(check_ott_id_size(i));
    }
    return ids;
}

vector<OttId> get_ids_from_file(const string& filename) {
    if (filename == "-") {
        return get_ids_from_stream(std::cin);
    }
    std::ifstream file(filename);
    return get_ids_from_stream(file);
}

vector<OttId> get_ids_from_tree(const string& filename) {
    vector<OttId> ids;
    auto tree = get_tree<Tree_t>(filename);
    for(auto nd: iter_post_const(*tree)) {
        auto id = nd->get_ott_id();
        if (id != -1) {
            ids.push_back(id);
        }
    }
    return ids;
}

vector<OttId> get_ids_matching_regex(const Taxonomy& taxonomy, const string& rgx) {
    std::regex e(rgx);
    vector<OttId> ids;
    for(const auto& rec: taxonomy) {
        std::cmatch m;
        if (std::regex_match(rec.name.data(), rec.name.data()+rec.name.size(), m, e)) {
            ids.push_back(rec.id);
        }
    }
    return ids;
}

bool has_flags(tax_flags flags, tax_flags any_flags, tax_flags all_flags) {
    if ((flags&all_flags) != all_flags) {
        return false;
     }
    if (any_flags.any() and (flags&any_flags).none()) {
        return false;
    }
    return true;
}

void flag_extinct_clade_as_incertae_sedis(Taxonomy & taxonomy) {
    auto nodeNamer = [](const auto& record){return string(record.name)+"_ott"+std::to_string(record.id);};
    auto tree = taxonomy.get_tree<RichTaxTree>(nodeNamer);
    const auto extinct_f = flags_from_string("extinct,extinct_inherited");
    set<OttId> not_extinct;
    set<OttId> to_add_incert;
    set<OttId> to_add_extinct;
    for (auto nd: iter_post(*tree)) {
        bool flag_par_not_extinct = false;
        const auto nd_id = nd->get_ott_id();
        if (contains(not_extinct, nd_id)) {
            flag_par_not_extinct = true;
        } else if ((nd->get_data().get_flags() & extinct_f).any()) {
            to_add_incert.insert(nd_id);
        } else {
            if (nd->get_first_child() == nullptr) {
                // tips not flagged as extinct are extant
                not_extinct.insert(nd_id);
                flag_par_not_extinct = true;
            } else {
                // internals will have been added to not_extinct if any child was extant
                to_add_incert.insert(nd_id);
                to_add_extinct.insert(nd_id);
            }
        }
        if (flag_par_not_extinct) {
            auto par_ptr = nd->get_parent();
            if (par_ptr == nullptr) {
                continue;
            }
            const auto par_id = par_ptr->get_ott_id();
            not_extinct.insert(par_id);
        }
    }
    

    const auto just_extinct_f = flags_from_string("extinct");
    const auto incert_sed_f = flags_from_string("incertae_sedis");
    for (auto& rec: taxonomy) {
        //std::cout << rec.name << "  orig = " << flags_to_string(rec.flags) << '\n';
        if (contains(to_add_incert, rec.id)) {
            if ((rec.flags & incert_sed_f).any()) {
                //std::cout << "already flagged inc-sed\n" ;
            } else {
                //std::cout << "adding inc-sed flag\n" ;
                rec.flags |= incert_sed_f;
                if (contains(to_add_extinct, rec.id)) {
                    std::cerr << "adding extinct flag to id=" << rec.id << "\n" ;
                    rec.flags |= just_extinct_f;
                }
            }
        } else {
            const auto eflags_present = rec.flags & extinct_f;
            if (eflags_present.any()) {
                std::cerr << "removing extinct flag from id=" << rec.id << "\n" ;
                rec.flags ^= eflags_present;
            } else {
                //std::cout << "not modifying flags\n" ;
            }
        }
        //std::cout << rec.name << "  final = " << flags_to_string(rec.flags) << '\n';
    }
}

void show_taxonomy_ids(const Taxonomy& taxonomy,
                       const string& format,
                       const vector<OttId>& ids,
                       std::function<bool(tax_flags)> flags_match) {
    for(auto id: ids) {
        try {
            auto& rec = taxonomy.record_from_id(id);
            if (flags_match(rec.flags)) {
                std::cout << format_with_taxonomy("No original label",format,rec,taxonomy) << "\n";
            }
        } catch (...) {
            std::cerr << "id=" << id << ": not in taxonomy!\n";
        }
    }
}

std::function<bool(tax_flags)> get_flags_match(variables_map& args) {
    tax_flags all_flags;
    if (args.count("all-flags")) {
        all_flags = flags_from_string(args["all-flags"].as<string>());
    }
    tax_flags any_flags;
    if (args.count("any-flags")) {
        any_flags = flags_from_string(args["any-flags"].as<string>());
    }
    if (any_flags.any() and all_flags.any()) {
        return [all_flags,any_flags](tax_flags flags) {return (flags&any_flags).any() and
                                                      (flags&all_flags)==all_flags; };
    } else if (any_flags.any()) {
        return [any_flags](tax_flags flags) { return (flags&any_flags).any(); };
    } else if (all_flags.any()) {
        return [all_flags](tax_flags flags) { return (flags&all_flags)==all_flags; };
    } else {
        if (args.count("cull-flags")) {
            tax_flags cull_flags;
            cull_flags = flags_from_string(args["cull-flags"].as<string>());
            return [cull_flags](tax_flags flags) { return ! ((flags&cull_flags).any()); };
        }
        return [](tax_flags){return true;};
    }
}

void report_dist_to_root(std::ostream & out, const Taxonomy & taxonomy) {
    for (auto & rec : taxonomy) {
        out << rec.id << '\t' << rec.depth << '\n';
    }
}

int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc, argv);
        auto format = args["format"].as<string>();
        auto flags_match = get_flags_match(args);
        auto taxonomy = load_taxonomy(args);
        const bool root_changed = args.count("root") || args.count("xroot");
        const bool detach_root = bool(args.count("xroot"));
        if (detach_root) {
            taxonomy[0].parent_id = 0;
        }
        if (args.count("extinct-to-incert")) {
            flag_extinct_clade_as_incertae_sedis(taxonomy);
        }
        if (args.count("show-root")) {
            show_rec(taxonomy[0]);
            return 0;
        } else if (args.count("in-tree")) {
            auto ids = get_ids_from_tree(args["in-tree"].as<string>());
            show_taxonomy_ids(taxonomy, format, ids, flags_match);
            return 0;
        } else if (args.count("in-file")) {
            auto ids = get_ids_from_file(args["in-file"].as<string>());
            show_taxonomy_ids(taxonomy, format, ids, flags_match);
            return 0;
        } else if (args.count("find")) {
            vector<OttId> ids = get_ids_matching_regex(taxonomy, args["find"].as<string>());
            show_taxonomy_ids(taxonomy, format, ids, flags_match);
            return 0;
        } else if (args.count("id")) {
            OttId id = args["name"].as<OttId>();
            show_taxonomy_ids(taxonomy, format, {id}, flags_match);
            return 0;
        } else if (args.count("any-flags") or args.count("all-flags")) {
            string format=args["format"].as<string>();
            for(const auto& rec: taxonomy) {
                if (flags_match(rec.flags)) {
                    std::cout << format_with_taxonomy("No original label",format,rec,taxonomy) << "\n";
                }
            }
            return 0;
        } else if (args.count("cull-flags")) {
            string format=args["format"].as<string>();
            for(const auto& rec: taxonomy) {
                if (flags_match(rec.flags)) {
                    std::cout << format_with_taxonomy("No original label",format,rec,taxonomy) << "\n";
                }
            }
            return 0;
        }
        if (args.count("degree")) {
            long id = args["degree"].as<long>();
            std::cout << "degree = " << taxonomy.record_from_id(id).out_degree << std::endl;
            return 0;
        } else if (args.count("children")) {
            long id = args["children"].as<long>();
            for(const auto& rec: taxonomy) {
                OttId parent_id = taxonomy[rec.parent_index].id;
                if (parent_id == id) {
                    show_rec(rec);
                }
            }
            return 0;
        } else if (args.count("parent")) {
            OttId id = args["parent"].as<OttId>();
            auto parent_index = taxonomy.record_from_id(id).parent_index;
            show_rec(taxonomy[parent_index]);
            return 0;
        }
        else if (args.count("high-degree-nodes")) {
            auto n = args["high-degree-nodes"].as<int>();
            auto index_vec = get_index_vec(taxonomy.size());
            std::sort(index_vec.begin(),
                      index_vec.end(),
                      [&taxonomy](int i, int j) {
                        return taxonomy[i].out_degree > taxonomy[j].out_degree;
                      });
            for(long i = 0; i < n; i++) {
                auto id = taxonomy[index_vec[i]].id;
                show_rec(taxonomy.record_from_unforwarded_id(id));
            }
            return 0;
        } else if (args.count("write-tree")) {
            auto nodeNamer = [](const auto& record){return string(record.name)+"_ott"+std::to_string(record.id);};
            write_tree_as_newick(cout, *taxonomy.get_tree<Tree_t>(nodeNamer));
            std::cout << std::endl;
        }
        if (args.count("write-taxonomy")) {
            taxonomy.write(args["write-taxonomy"].as<string>(), false, !root_changed);
        }
        if (args.count("name")) {
            OttId id = args["name"].as<OttId>();
            std::cout << taxonomy.record_from_id(id).name << std::endl;
        }
        if (args.count("uniqname")) {
            OttId id = args["uniqname"].as<OttId>();
            std::cout << taxonomy.record_from_id(id).uniqname << std::endl;
        }
        if (args.count("report-lost-taxa")) {
            string treefile = args["report-lost-taxa"].as<string>();
            report_lost_taxa(taxonomy,treefile);
        }
        if (args.count("version")) {
            std::cout << taxonomy.get_version() << std::endl;
        }
        if (args.count("report-dist-to-root")) {
            report_dist_to_root(std::cout, taxonomy);
        }
    } catch (std::exception& e) {
        cerr << "otc-taxonomy-parser: Error! " << e.what() << std::endl;
        return 1;
    }
}

// 1. Write a parser to read the lines faster
// 2. Avoid memory allocation -- by mmapping the taxonomy file?
// 3. Convert the flags into a bitmask
// 4. Should the Rank be a converted to an integer?
// 5. Can we assign OTT IDs to internal nodes of a tree while accounting for Incertae Sedis taxa?
// * What are the triplet-inference rules for the Incertae Sedis problem?

// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

