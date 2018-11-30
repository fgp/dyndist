/*
 * dyndist/vector_distribution.h
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

#ifndef dyndist_vector_distribution_h
#define dyndist_vector_distribution_h

#include <utility>
#include <iterator>
#include <vector>
#include <algorithm>

#include "dyndist/global.h"
#include "dyndist/discrete_distribution.h"

DYNDIST_NAMESPACE_BEGIN

namespace vector_distribution_details
{
    template<typename IteratorType, typename EventsIteratorType>
    struct iterator_base {
    public:
        typedef IteratorType iterator_type;

        typedef EventsIteratorType event_iterator_type;

    private:
        iterator_type* self()
        { return static_cast<iterator_type*>(this); }

        const iterator_type* self() const
        { return static_cast<const iterator_type*>(this); }

    public:
        explicit iterator_base(const event_iterator_type event_i)
            :m_event_i(event_i)
        {}

        iterator_type& operator++()
        { ++m_event_i; return *self(); }

        iterator_type& operator+=(std::ptrdiff_t d)
        { m_event_i += d; return *self(); }

        iterator_type& operator--()
        { --m_event_i; return *self(); }

        iterator_type& operator-=(std::ptrdiff_t d)
        { m_event_i -= d; return *self(); }

        iterator_type operator++(int)
        { iterator_type i = *self(); ++(*self()); return i; }

        iterator_type operator+(std::ptrdiff_t d) const
        { iterator_type i = *self(); i += d; return i; }

        iterator_type operator--(int)
        { iterator_type i = *self(); --(*self()); return i; }

        iterator_type operator-(std::ptrdiff_t d) const
        { iterator_type i = *self(); i -= d; return i; }

        std::ptrdiff_t operator-(const iterator_type& o) const
        { return m_event_i - o.m_event_i; }

        bool operator< (const iterator_type& o) const
        { return m_event_i < o.m_event_i; }

        bool operator<= (const iterator_type& o) const
        { return m_event_i <= o.m_event_i; }

        bool operator==(const iterator_type& o) const
        { return m_event_i == o.m_event_i; }

        bool operator!=(const iterator_type& o) const
        { return m_event_i != o.m_event_i; }

        bool operator>=(const iterator_type& o) const
        { return m_event_i >= o.m_event_i; }

        bool operator> (const iterator_type& o) const
        { return m_event_i > o.m_event_i; }

    protected:
        event_iterator_type m_event_i;
    };
}

/**
 * \class vector_distribution
 *
 * Implements an distribution over events indexed by 0,...,n-1 with arbitrary weights
 *
 * \tparam WeightType The integral type used for represent element weights
 */
template<typename WeightType>
class vector_distribution
{
public:
    typedef WeightType weight_t;

private:
    typedef discrete_distribution_pointer<std::size_t, weight_t, std::size_t>
            event_pointer;

    typedef std::vector<event_pointer>
            events_vector_type;

    /** The EventMovedFunctor for distribution_type */
    struct species_moved
    {
        explicit species_moved(events_vector_type& vector)
            :m_vector(vector)
        {}

        DYNDIST_INLINE
        void operator()(event_pointer species, event_pointer species_new)
        {
            DYNDIST_ASSERT(species->data == species_new->data);
            const std::size_t k = species->data;
            DYNDIST_ASSERT(k < m_vector.size());
            DYNDIST_ASSERT(m_vector[k] == species);
            m_vector[k] = species_new;
        }

    private:
        events_vector_type& m_vector;
    };

public:
    typedef discrete_distribution<std::size_t, weight_t, std::size_t, species_moved>
            distribution_type;

    struct iterator;

private:
    typedef typename distribution_type::value_type event_type;

    /** Proxy returned by operator= to allow new weights to be assigned */
    struct event_proxy
    {
        /** attempting to take the address of an event_proxy returns an iterator pointing to the event */
        DYNDIST_INLINE_FLATTEN
        iterator operator& () const
        { return iterator(m_vector_distribution, m_event_idx); }

        DYNDIST_INLINE_FLATTEN
        operator weight_t () const
        { return m_vector_distribution.at(m_event_idx); }

        template<typename T>
        DYNDIST_INLINE
        event_proxy& operator=(const T& w)
        { m_vector_distribution.replace(m_event_idx, w); return *this; }

        template<typename T>
        DYNDIST_INLINE
        const T& operator+=(const T& w)
        { m_vector_distribution.replace(m_event_idx, m_vector_distribution.at(m_event_idx) + w); return w; }

        template<typename T>
        DYNDIST_INLINE
        const T& operator-=(const T& w)
        { m_vector_distribution.replace(m_event_idx, m_vector_distribution.at(m_event_idx) - w); return w; }

        template<typename T>
        DYNDIST_INLINE
        const T& operator*=(const T& w)
        { m_vector_distribution.replace(m_event_idx, m_vector_distribution.at(m_event_idx) * w); return w; }

        template<typename T>
        DYNDIST_INLINE
        const T& operator/=(const T& w)
        { m_vector_distribution.replace(m_event_idx, m_vector_distribution.at(m_event_idx) / w); return w; }

    private:
        friend class vector_distribution;
        friend class iterator;

        DYNDIST_INLINE
        event_proxy(vector_distribution& idx_dist, std::size_t event_idx)
                :m_vector_distribution(idx_dist)
                ,m_event_idx(event_idx)
        {}

        vector_distribution& m_vector_distribution;
        const std::size_t m_event_idx;
    };

public:
    struct const_iterator;

    /** Mutable iterator for vector_distribution */
    struct iterator : vector_distribution_details::iterator_base<iterator, std::size_t>
    {
        DYNDIST_INLINE_FLATTEN
        event_proxy operator*() const
        { return event_proxy(*m_vector_distribution, this->m_event_i); }

        DYNDIST_INLINE_FLATTEN
        const weight_t* operator->() const
        { return &((*m_vector_distribution)[this->m_event_i]); }

    private:
        friend class const_iterator;
        friend class vector_distribution;

        DYNDIST_INLINE
        iterator(vector_distribution& distribution, std::size_t event_idx)
                :vector_distribution_details::iterator_base<iterator, std::size_t>(event_idx)
                ,m_vector_distribution(&distribution)
        {}

        vector_distribution* m_vector_distribution;
    };

    /** Const iterator for vector_distribution */
    struct const_iterator : vector_distribution_details::iterator_base<const_iterator, typename events_vector_type::const_iterator>
    {
        /** Creates a const_iterator from a non-const iterator */
        DYNDIST_INLINE_FLATTEN
        const_iterator(const iterator& it)
            :vector_distribution_details::iterator_base<const_iterator, typename events_vector_type::const_iterator>(
                    it.m_vector_distribution->m_events.begin() + it.m_event_i)
        {}

        DYNDIST_INLINE_FLATTEN
        const weight_t& operator*() const
        { return (*this->m_event_i)->weight; }

        DYNDIST_INLINE_FLATTEN
        const weight_t* operator->() const
        { return &(*this->m_event_i)->weight; }

    private:
        friend class vector_distribution;

        DYNDIST_INLINE_FLATTEN
        explicit const_iterator(typename events_vector_type::const_iterator event_i)
                :vector_distribution_details::iterator_base<const_iterator, typename events_vector_type::const_iterator>(event_i)
        {}
    };

public:
    /** Constructs an empty distribution */
    vector_distribution()
        :m_distribution(species_moved(m_events))
    {}

    /** Constructs a distribution containing \param count elements of weight \param weight */
    vector_distribution(std::size_t count, std::size_t weight);

    /** Returns the number of events in the distribution */
    DYNDIST_INLINE
    std::size_t size() const
    { return m_events.size(); }

    /** Returns the total weight of the events in the distribution */
    DYNDIST_INLINE
    weight_t weight() const
    { return m_distribution.weight(); }

    /** Returns the weight of the i-th event */
    DYNDIST_INLINE
    const weight_t& at(std::size_t i) const
    { return m_events.at(i)->weight; }

    /** Updates the weight of the i-th event */
    weight_t replace(std::size_t i, const weight_t& w)
    { return m_distribution.set(m_events.at(i), w); }

    /** Returns the weight of the i-th event */
    DYNDIST_INLINE
    const weight_t& operator[](std::size_t i) const
    { return at(i); }

    /** Returns a mutable proxy of the i-th event that allows a new weight to be assigned with = */
    DYNDIST_INLINE
    event_proxy operator[](std::size_t i)
    { return event_proxy(*this, i); }

    /** Adds a event with weight w, returns the event's index */
    std::size_t push_back(const weight_t& w);

    /** Draws a random event */
    template<typename Engine>
    std::size_t operator()(Engine& engine);

    DYNDIST_INLINE
    iterator begin()
    { return iterator(*this, 0); }

    DYNDIST_INLINE
    const_iterator begin() const
    { return cbegin(); }

    DYNDIST_INLINE
    const_iterator cbegin()
    { return const_iterator(m_events.cbegin()); }

    DYNDIST_INLINE
    iterator end()
    { return iterator(*this, m_events.size()); }

    DYNDIST_INLINE
    const_iterator end() const
    { return cend(); }

    DYNDIST_INLINE
    const_iterator cend()
    { return const_iterator(m_events.cend()); }

private:
    events_vector_type m_events;
    distribution_type m_distribution;
};

template <typename W>
vector_distribution<W>::vector_distribution(std::size_t count, std::size_t weight)
    :m_distribution(species_moved(m_events))
{
    for(std::size_t i = 0; i < count; ++i)
        push_back(weight);
}

/** Adds a event with weight w, returns the event's index */
template<typename W>
std::size_t
DYNDIST_FLATTEN
vector_distribution<W>::push_back(const weight_t& w)
{ const std::size_t i = m_events.size(); m_events.push_back(m_distribution.insert(event_type(w, i))); return i;}

/** Draws a random event */
template<typename W>
template<typename Engine>
DYNDIST_FLATTEN
std::size_t
vector_distribution<W>::operator()(Engine& engine)
{ return m_distribution(engine)->data; }

DYNDIST_NAMESPACE_END

#endif
