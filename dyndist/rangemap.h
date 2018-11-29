/*
 * dyndist/rangemap.h
 *
 * Copyright (c) 2014,2018 Florian Pflug.
 *
 * This file is part of dyndist
 *
 * dyndist is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * dyndist is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with dyndist.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef dyndist_rangemap_h
#define dyndist_rangemap_h

#include "dyndist/global.h"

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <deque>

DYNDIST_NAMESPACE_BEGIN

template<typename ValueType>
class rangemap : std::deque< std::pair<std::ptrdiff_t, ValueType> >
{
public:
    typedef std::ptrdiff_t key_type;
    typedef ValueType mapped_type;

private:
    typedef std::deque< std::pair<std::ptrdiff_t, ValueType> > super_type;
    
    /* Hide non-applicable members */
    void resize();
    void erase();
    
public:
    using super_type::size;
    using super_type::front;
    using super_type::back;
    using super_type::empty;
    using super_type::begin;
    using super_type::end;
    using super_type::rbegin;
    using super_type::rend;
    using super_type::pop_front;
    using super_type::pop_back;
    using super_type::clear;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::iterator iterator;
    typedef typename super_type::const_iterator const_iterator;
    typedef typename super_type::reverse_iterator reverse_iterator;
    typedef typename super_type::const_reverse_iterator const_reverse_iterator;
    
    key_type lower_bound() const {
        DYNDIST_ASSERT(!empty());
        return front().first;
    }
    
    key_type upper_bound() const {
        DYNDIST_ASSERT(!empty());
        return back().first;
    }
    
    bool contains(const key_type key) const {
        return !empty() && (lower_bound() <= key) && (key <= upper_bound());
    }
    
    mapped_type& operator[](const key_type key) {
        DYNDIST_ASSERT(contains(key));
        return super_type::operator[]((std::size_t)(key - lower_bound())).second;
    }
    const mapped_type& operator[](const key_type key) const {
        DYNDIST_ASSERT(contains(key));
        return super_type::operator[](key - lower_bound()).second;
    }

    mapped_type& at(const key_type key) {
        if (!contains(key))
            throw std::out_of_range("no such key");
        return super_type::operator[](key - lower_bound()).second;
    }
    const mapped_type& at(const key_type key) const {
        if (!contains(key))
            throw std::out_of_range("no such key");
        return super_type::operator[](key - lower_bound()).second;
    }
    
    iterator locate(const key_type key) {
        DYNDIST_ASSERT(contains(key));
        DYNDIST_ASSERT((begin() + (key - lower_bound()))->first == key);
        return begin() + (key - lower_bound());
    }
    const_iterator locate(const key_type key) const {
        DYNDIST_ASSERT(contains(key));
        DYNDIST_ASSERT((begin() + (key - lower_bound()))->first == key);
        return begin() + (key - lower_bound());
    }
    
    iterator find(key_type key);
    const_iterator find(key_type key) const;
    
    value_type& push_front(const mapped_type& mapped_value = mapped_type());

    value_type& push_back(const mapped_type& mapped_value = mapped_type());

    template<typename F>
    void enlarge_to(key_type lower_bound, key_type upper_bound, F& callback);

    void restrict_to(key_type upper_bound, key_type lower_bound);

    std::pair<iterator,bool> insert(const value_type& value);
    
    iterator insert(iterator /* hint */, const value_type& value) {
        /* Hint is ignored */
        return insert(value).second;
    }
    
    template<class InputIterator>
    void insert(InputIterator first, InputIterator last) {
        for(; first != last; ++first)
            insert(*first);
    }
};

template<typename V>
inline
typename rangemap<V>::iterator
rangemap<V>::find(const key_type key)
{
    if (empty())
        return end();
    
    const std::ptrdiff_t d_lb = key - lower_bound();
    if ((d_lb < 0) || (upper_bound() < key))
        return end();
    else {
        DYNDIST_ASSERT((begin() + d_lb)->first == key);
        return begin() + d_lb;
    }
}
template<typename V>
inline
typename rangemap<V>::const_iterator
rangemap<V>::find(const key_type key) const
{
    if (empty())
        return end();
    
    const std::ptrdiff_t d_lb = key - lower_bound();
    if ((d_lb < 0) || (upper_bound() < key))
        return end();
    else {
        DYNDIST_ASSERT((begin() + d_lb)->first == key);
        return begin() + d_lb;
    }
}

template<typename V>
typename rangemap<V>::value_type&
rangemap<V>::push_front(const mapped_type& mapped_value)
{
    DYNDIST_ASSERT(!empty());
    super_type::push_front(std::make_pair(lower_bound()-1, mapped_value));
    return front();
}

template<typename V>
typename rangemap<V>::value_type&
rangemap<V>::push_back(const mapped_type& mapped_value)
{
    DYNDIST_ASSERT(!empty());
    super_type::push_back(std::make_pair(upper_bound()+1, mapped_value));
    return back();
}

template<typename V> template<typename F>
void
rangemap<V>::enlarge_to(const key_type lb, const key_type ub, F& callback)
{
    if (empty()) {
        super_type::push_front(std::make_pair(lb, mapped_type()));
        callback(front(), true);
    }

    key_type curr_lb = lower_bound();
    for(; lb < curr_lb; --curr_lb) {
        super_type::push_front(std::make_pair(curr_lb - 1, mapped_type()));
        callback(front(), true);
    }
    
    key_type curr_ub = upper_bound();
    for(; ub > curr_ub; ++curr_ub) {
        super_type::push_back(std::make_pair(curr_ub + 1, mapped_type()));
        callback(back(), false);
    }
    
    DYNDIST_ASSERT((lower_bound() <= lb) && (ub <= upper_bound()));
}

template<typename V>
void
rangemap<V>::restrict_to(const key_type lb, const key_type ub)
{
    if (empty())
        return;
    
    key_type curr_lb = lower_bound(), curr_ub = upper_bound();
    DYNDIST_ASSERT(lb <= ub);
    if ((ub < curr_lb) || (curr_ub < lb)) {
        clear();
        return;
    }
    if (curr_lb < lb)
        erase(begin(), begin() + (lb - curr_lb));
    if (ub < curr_ub)
        erase(end() - (curr_ub - ub), end());
    
    DYNDIST_ASSERT((lb <= lower_bound()) && (upper_bound() <= ub));
}

template<typename V>
std::pair< typename rangemap<V>::iterator, bool >
rangemap<V>::insert(const value_type& value)
{
    const key_type& key = value.first;
    if (empty()) {
        super_type::push_front(value);
        return std::make_pair(begin(), true);
    }

    key_type lb = lower_bound(), ub = upper_bound();
    for(; key < lb; --lb) {
        if (key == lb - 1) {
            super_type::push_front(value);
            return std::make_pair(begin(), true);
        }
        super_type::push_front(std::make_pair(lb - 1, mapped_type()));
    }
    for(; key > ub; ++ub) {
        if (key == ub + 1) {
            super_type::push_back(value);
            return std::make_pair(end() - 1, true);
        }
        super_type::push_back(std::make_pair(ub + 1, mapped_type()));
    }

    DYNDIST_ASSERT(contains(key));
    DYNDIST_ASSERT((begin() + (key - lb))->first == (end() - (ub - key) - 1)->first);
    return std::make_pair(begin() + (key - lb), false);
}

DYNDIST_NAMESPACE_END

#endif
