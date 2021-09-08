#ifndef CACHE_LIST
#define CACHE_LIST

#include <list>
#include "otc/assert.hh"
#include <initializer_list>
#include <boost/container/list.hpp>
// Currently, since I'm not using T*, this only works for POD type.
// Perhaps I should clear list entries by setting their x equal to T()?

template <typename T>
class cached_list
{
    using base_list = boost::container::list<T>;

    /// This points to the first free node.
    static base_list free_pool;

    /// The actual storage.
    base_list L;

public:
    typedef typename base_list::iterator iterator;
    typedef typename base_list::const_iterator const_iterator;

    int size() const {return L.size();}

    const_iterator begin() const 
        {
            return L.begin();
        }

    iterator begin()
        {
            return L.begin();
        }

    const_iterator end() const 
        {
            return L.end();
        }

    iterator end() {
        return L.end();
    }

    /// Remove an element from the list. (Place it in the free area.)
    void erase(iterator loc) 
        {
            // Move (used) element from loc to the free pool.
            free_pool.splice(free_pool.begin(), L, loc);
        }

    /// Insert the value t *before* loc.
    iterator insert(const_iterator loc, const T& t)
        {
            if (free_pool.empty())
                return L.insert(loc, t);
            else
            {
                iterator free_loc = free_pool.begin();
      
                // Set the value
                *free_loc = t;

                // Move (used) element from loc to *before* begin_free.
                L.splice(loc, free_pool, free_loc);

                return free_loc;
            }
        }

    void splice(const_iterator loc, cached_list<T>& CL2)
        {
            L.splice(loc, CL2.L);
        }

    void splice(const_iterator loc, cached_list<T>& CL2, const_iterator loc2a, const_iterator loc2b)
        {
            L.splice(loc, CL2.L, loc2a, loc2b);
        }

    const T& back() const
        {
            return L.back();
        }

    T& back()
        {
            return L.back();
        }

    void pop_back()
        {
            assert(not empty());

            iterator last = end();
            last--;
            erase(last);
        }

    iterator push_back(const T& t)
        {
            return insert(end(), t);
        }

    void pop_front()
        {
            assert(not empty());

            erase(begin());
        }

    iterator push_front(const T& t)
        {
            return insert(begin(), t);
        }

    bool empty() const
        {
            return L.empty();
        }

    int count(const T& t) const
        {
            int total = 0;
            for(const auto& i: L)
                if (i == t)
                    total++;
            return total;
        }

    void clear() 
        {
            free_pool.splice(free_pool.begin(), L);
        }

    void swap(cached_list<T>& CL2) noexcept
        {
            L.swap(CL2.L);
        }

    cached_list<T>& operator=(const base_list& L2)
        {
            clear();
            for(auto i: L2)
                push_back(i);
            return *this;
        }

    cached_list<T>& operator=(const cached_list<T>& CL2)
        {
            return operator=(CL2.L);
        }

    cached_list<T>& operator=(base_list&& L2) noexcept
        {
            clear();
            L.swap(L2);
            return *this;
        }

    cached_list<T>& operator=(cached_list<T>&& CL2) noexcept
        {
            return operator=(std::move(CL2.L));
        }

    cached_list() {}

    cached_list(const base_list& L2)
        {
            operator=(L2);
        }

    cached_list(const base_list&& L2)
        { 
            operator=(std::move(L2));
        }

    cached_list(std::initializer_list<T> IL2)
        { 
            for(auto i: IL2)
                push_back(i);
        }

    cached_list(const cached_list<T>& L)
        :cached_list()
        {
            operator=(L);
        }

    cached_list(cached_list<T>&& L) noexcept
        :cached_list()
        { 
            swap(L);
        }

    ~cached_list()
        {
            clear();
        }
};

template<typename T> typename cached_list<T>::base_list cached_list<T>::free_pool = cached_list<T>::base_list{};

#endif
