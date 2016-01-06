#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

struct taxonomy_record
{
    int id = 0;
    int parent_id = 0;
    string name;
    taxonomy_record(int i1, int i2, string&& s): id(i1), parent_id(i2), name(std::move(s)) {}
    taxonomy_record(int i1, int i2, const string& s): id(i1), parent_id(i2), name(s) {}
};

int main(int argc, char* argv[])
{
    try
    {
        if (argc != 2)
            throw OTCError()<<"Expecting exactly 1 argument, but got "<<argc-1<<".";

        string filename = argv[1];

        std::ifstream taxonomy(filename);
        if (not taxonomy)
            throw OTCError()<<"Could not open file '"<<filename<<"'.";
    
        string line;
        int count = 0;
        vector<taxonomy_record> lines;
        std::unordered_map<long,Tree_t::node_type*> id_to_node;
        std::getline(taxonomy,line);
        if (line != "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t")
            throw OTCError()<<"First line of file '"<<filename<<"' is not a taxonomy header.";

        std::unique_ptr<Tree_t> tree(new Tree_t());
        while(std::getline(taxonomy,line))
        {
            int b1 = line.find("\t|\t",0);
            const char* start = line.c_str();
            char* end = nullptr;
            int id = std::strtoul(line.c_str(),&end,10);
            start = end + 3;
            int parent_id = std::strtoul(start,&end,10);
//            lines.push_back(line);
            start = end + 3;
            const char* end3 = std::strstr(start,"\t|\t");
            string name = line.substr(start - line.c_str(),end3 - start);
//            cerr<<id<<"\t"<<parent_id<<"\t'"<<name<<"'\n";

            Tree_t::node_type* nd = nullptr;
            if (not parent_id)
            {
                nd = tree->createRoot();
            }
            else
            {
                auto parent_nd = id_to_node.at(parent_id);
                nd = tree->createChild(parent_nd);
            }
            nd->setOttId(id);
            nd->setName(name);
            id_to_node[id] = nd;
            
            lines.push_back({id,parent_id,name});
            count++;
        }
        cerr<<"#lines = "<<count<<std::endl;
        //      writeTreeAsNewick(cout, *tree);
    }
    catch (std::exception& e)
    {
        cerr<<"otc-taxonomy-parser: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}
