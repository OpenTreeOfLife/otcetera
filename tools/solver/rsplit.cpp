#include "tools/solver/rsplit.h"
#include "otc/util.h"

using std::set;
using std::vector;

RSplitObj::RSplitObj()
{
    id = num++;
}

RSplitObj::RSplitObj(const set<int>& i, const set<int>& a)
{
    id = num++;
    in  = set_to_vector(i);
    set_difference(begin(a), end(a), begin(in), end(in), std::inserter(out, out.end()));
}

RSplitObj::RSplitObj(const set<int>& i, const vector<int>& a)
{
    id = num++;
    in  = set_to_vector(i);
    set_difference(begin(a), end(a), begin(in), end(in), std::inserter(out, out.end()));
}

std::size_t RSplitObj::num = 0;

std::ostream& operator<<(std::ostream& o, const RSplitObj& s)
{
    o<<"["<<s.in.size() + s.out.size()<<" tips] ";
    for(auto& is: s.in)
        o<<"ott"<<is<<" ";
    o<<"| ";
    for(auto& os: s.out)
        o<<"ott"<<os<<" ";
    return o;
}

RSplit split_from_include_exclude(const set<int>& i, const set<int>& e)
{
    RSplit s(new RSplitObj);
    s->in = set_to_vector(i);
    s->out = set_to_vector(e);
    return s;
}

std::ostream& operator<<(std::ostream& o, const ConstRSplit& s)
{
    otc::write_separated_collection(o, s->in, " ") <<" | ";
    if (s->out.size() < 100)
        otc::write_separated_collection(o, s->out, " ");
    else
    {
        auto it = s->out.begin();
        for(int i=0;i<100;i++)
            o << *it++ <<" ";
        o << "...";
    }
    return o;
}

std::ostream& operator<<(std::ostream& o, const RSplit& s)
{
    return o<<(ConstRSplit(s));
}
