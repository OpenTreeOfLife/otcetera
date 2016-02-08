// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <bitset>
#include <fstream>

#include <boost/filesystem/operations.hpp>

namespace fs = boost::filesystem;

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::ofstream;

using boost::string_ref;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace otc
{

    taxonomy_record::taxonomy_record(const string& line_)
        :line(line_)
    {
        // parse the line
        // also see boost::make_split_iterator
        const char* start[8];
        const char* end[8];

        start[0] = line.c_str();
        for(int i=0;i<7;i++)
        {
            end[i] = std::strstr(start[i],"\t|\t");
            start[i+1] = end[i] + 3;
        }

//    boost::container::small_vector<string_ref,10> words;
//    for (string_ref&& r : iter_split(v, line, token_finder(is_any_of(","))) |
//             transformed([](R const& r){return boost::string_ref(&*r.begin(), r.size());})
//        )
//        words.push_back(r);

        char *temp;
        id = std::strtoul(start[0], &temp, 10);
        parent_id = std::strtoul(start[1], &temp, 10);
        name = string_ref(start[2],end[2]-start[2]);
        rank = string_ref(start[3],end[3]-start[3]);
        sourceinfo = string_ref(start[4],end[4]-start[4]);
        uniqname = string_ref(start[5],end[5]-start[5]);
        flags = flags_from_string(start[6],end[6]);
    }

    long Taxonomy::map(long old_id) const
    {
        if (index.count(old_id))
            return old_id;
            
        auto loc = forwards.find(old_id);
        if (loc != forwards.end())
            return loc->second;

        return -1;
    }
    
    void Taxonomy::write(const std::string& newdirname)
    {
        fs::path old_dir = path;
        fs::path new_dir = newdirname;

        if (fs::exists(new_dir))
            throw OTCError()<<"File '"<<newdirname<<"' already exists!";

        fs::create_directories(new_dir);

        // Copy the other files.
        for(const auto& name: {"about.json", "conflicts.tsv", "deprecated.tsv",
                    "log.tsv", "otu_differences.tsv", "synonyms.tsv", "weaklog.csv"})
            fs::copy_file(old_dir/name,new_dir/name);

        // Write the new version file.
        {
            ofstream version_file((new_dir/"version.txt").string());
            version_file<<version;

            if (keep_root != -1 or cleaning_flags.any())
            {
                version_file<<"modified: ";
                if (keep_root != -1)
                    version_file<<"root="<<keep_root<<"  ";
                if (cleaning_flags.any())
                    version_file<<"cleaning_flags="<<flags_to_string(cleaning_flags);
                version_file<<"\n";
            }
            version_file.close();
        }

        // Write the new taxonomy file.
        {
            ofstream tf ((new_dir/"taxonomy.tsv").string());
            tf << "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t"<<std::endl;
            string sep = "\t|\t";
            for(const auto& r: *this)
                tf << r.line <<"\n";
            tf.close();
        }

        // Write the new forwards file.
        {
            ofstream ff((new_dir/"forwards.tsv").string());
            ff << "id\treplacement\n";
            for(const auto& p: forwards)
                ff<<p.first<<'\t'<<p.second<<'\n';
            ff.close();
        }
    }
    
    Taxonomy::Taxonomy(const string& dir, bitset<32> cf, int kr)
        :keep_root(kr),
        cleaning_flags(cf),
         path(dir),
         version(readStrContentOfUTF8File(dir + "/version.txt"))
    {
        string filename = path + "/taxonomy.tsv";

        // 1. Open the file.
        std::ifstream taxonomy_stream(filename);
        if (not taxonomy_stream)
            throw OTCError()<<"Could not open file '"<<filename<<"'.";

        // 2. Read and check the first line
        string line;
        std::getline(taxonomy_stream,line);
        if (line != "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t")
            throw OTCError()<<"First line of file '"<<filename<<"' is not a taxonomy header.";

        // 3. Read records up to the record containing the root.
        int count = 0;
        if (keep_root != -1)
        {
            while(std::getline(taxonomy_stream,line))
            {
                count++;

                // Add line to vector
                emplace_back(line);

                if (back().id == keep_root)
                    break;

                pop_back();
            }
            if (empty())
                throw OTCError()<<"Root id '"<<keep_root<<"' not found.";
        }
        else
        {
            std::getline(taxonomy_stream,line);
            count++;
            // Add line to vector
            emplace_back(line);
        }
        if ((back().flags & cleaning_flags).any())
            throw OTCError()<<"Root taxon (ID = "<<back().id<<") removed according to cleaning flags!";
        index[back().id] = size() - 1;

        // 4. Read the remaining records
        while(std::getline(taxonomy_stream,line))
        {
            count++;

            // Add line to vector
            emplace_back(line);

            // Eliminate records that match the cleaning flags
            if ((back().flags & cleaning_flags).any())
            {
                pop_back();
                continue;
            }
        
            // Eliminate records whose parents have been eliminated, or are not found.
            auto loc = index.find(back().parent_id);
            if (loc == index.end())
            {
                pop_back();
                continue;
            }

            back().parent_index = loc->second;
        
            index[back().id] = size() - 1;
        }
        cerr<<"#lines read = "<<count<<std::endl;
        cerr<<"size = "<<size()<<std::endl;
        taxonomy_stream.close();

        std::ifstream forwards_stream(path + "/forwards.tsv");
        if (not forwards_stream)
            throw OTCError()<<"Could not open file '"<<filename<<"'.";

        std::getline(forwards_stream, line);
        while(std::getline(forwards_stream, line))
        {
            char* temp;
            long old_id = std::strtoul(line.c_str(), &temp, 10);
            assert(*temp == '\t');
            const char* temp2 = temp+1;
            long new_id = std::strtoul(temp2, &temp, 10);
            if (index.count(new_id))
                forwards[old_id] = new_id;
        }
    }

}
