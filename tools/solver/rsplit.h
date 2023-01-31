#ifndef RSPLIT_H
#define RSPLIT_H

#include <vector>
#include <optional>
#include <set>
#include <boost/smart_ptr/intrusive_ptr.hpp>

struct RSplitObj
{
    static std::size_t num;

    mutable int _refs = 0;

    friend inline void intrusive_ptr_release(RSplitObj* pThis)
    {
        if (--pThis->_refs == 0 ) {
            delete pThis;
        }
    }

    friend inline void intrusive_ptr_add_ref(RSplitObj* pThis)
    {
        pThis->_refs++;
    }

    friend inline void intrusive_ptr_release(const RSplitObj* pThis)
    {
        if(--const_cast<RSplitObj*>(pThis)->_refs == 0 ) {
            delete const_cast<RSplitObj*>(pThis);
        }
    }

    friend inline void intrusive_ptr_add_ref(const RSplitObj* pThis)
    {
        const_cast<RSplitObj*>(pThis)->_refs++;
    }

    std::vector<int> in;
    std::vector<int> out;
    std::optional<std::size_t> id;

    RSplitObj();
    RSplitObj(const std::set<int>& i, const std::set<int>& a);
    RSplitObj(const std::set<int>& i, const std::vector<int>& a);
};

// Reference-counted objects, with reference stored inside the object
using RSplit = boost::intrusive_ptr<RSplitObj>;
using ConstRSplit = boost::intrusive_ptr<const RSplitObj>;

std::ostream& operator<<(std::ostream& o, const RSplitObj& s);

std::ostream& operator<<(std::ostream& o, const ConstRSplit& s);

std::ostream& operator<<(std::ostream& o, const RSplit& s);

RSplit split_from_include_exclude(const std::set<int>& i, const std::set<int>& e);

/// Create a SORTED vector from a set
/// Better: myset | ranges::to<vector>
template <typename T>
std::vector<T> set_to_vector(const std::set<T>& s) {
    std::vector<T> v;
    v.reserve(s.size());
    std::copy(s.begin(), s.end(), std::back_inserter(v));
    return v;
}

#endif /* RSPLIT_H */
