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

private:
    typedef typename distribution_type::value_type event_type;

    /** Proxy returned by operator= to allow new weights to be assigned */
    struct event_proxy
    {
        event_proxy(vector_distribution& idx_dist, std::size_t event_idx)
            :m_vector_distribution(idx_dist)
            ,m_event_idx(event_idx)
        {}

        operator weight_t () const
        { return m_vector_distribution.at(m_event_idx); }

        template<typename T>
        event_proxy& operator=(const T& w)
        { m_vector_distribution.replace(m_event_idx, w); return *this; }

        template<typename T>
        const T& operator+=(const T& w)
        { m_vector_distribution.replace(m_event_idx, m_vector_distribution.at(m_event_idx) + w); return w; }

        template<typename T>
        const T& operator-=(const T& w)
        { m_vector_distribution.replace(m_event_idx, m_vector_distribution.at(m_event_idx) - w); return w; }

        template<typename T>
        const T& operator*=(const T& w)
        { m_vector_distribution.replace(m_event_idx, m_vector_distribution.at(m_event_idx) * w); return w; }

        template<typename T>
        const T& operator/=(const T& w)
        { m_vector_distribution.replace(m_event_idx, m_vector_distribution.at(m_event_idx) / w); return w; }

    private:
        vector_distribution& m_vector_distribution;
        const std::size_t m_event_idx;
    };

public:
    /** Constructs an empty distribution */
    vector_distribution()
        :m_distribution(species_moved(m_events))
    {}

    /** Constructs a distribution containing \param count elements of weight \param weight */
    vector_distribution(std::size_t count, std::size_t weight);

    /** Returns the number of events in the urn */
    std::size_t size() const
    { return m_events.size(); }

    weight_t weight() const
    { return m_distribution.weight(); }

    /** Returns the weight of the i-th event */
    const weight_t& at(std::size_t i) const
    { return m_events.at(i)->weight; }

    /** Updates the weight of the i-th event */
    weight_t replace(std::size_t i, const weight_t& w)
    { return m_distribution.set(m_events.at(i), w); }

    /** Returns the weight of the i-th event */
    const weight_t& operator[](std::size_t i) const
    { return at(i); }

    /** Returns a mutable proxy of the i-th event that allows a new weight to be assigned with = */
    event_proxy operator[](std::size_t i)
    { return event_proxy(*this, i); }

    /** Adds a event with weight w, returns the event's index */
    std::size_t push_back(const weight_t& w)
    { const std::size_t i = m_events.size(); m_events.push_back(m_distribution.insert(event_type(w, i))); return i;}

    /** Draws a random event */
    template<typename Engine>
    std::size_t operator()(Engine& engine)
    { return m_distribution(engine)->data; }

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


DYNDIST_NAMESPACE_END

#endif
