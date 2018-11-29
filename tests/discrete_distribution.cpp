/*
 * dyndist/tests/discrete_distribution.cpp
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

#define DYNDIST_ENABLE_VERIFICATION 1
#include "dyndist/global.h"

#include <cmath>
#include <limits>

#if __cplusplus < 201103L

/* C++ 98. Use boost for uint64_t and random number generation */
DYNDIST_NOWARN_PUSH
DYNDIST_NOWARN_SIGNS
#include <boost/integer.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
DYNDIST_NOWARN_POP
using boost::uint64_t;
using boost::random::mt19937;
using boost::random::uniform_int_distribution;
using boost::random::uniform_real_distribution;

#else // __cplusplus >= 201103L

/* C++ 11. Use standard library for uint64_t and random number generation */
#include <cstdint>
#include <random>
using std::uint64_t;
using std::mt19937;
using std::mt19937;
using std::uniform_int_distribution;
using std::uniform_real_distribution;

#endif // __cplusplus <> 201103L

DYNDIST_NOWARN_PUSH
DYNDIST_NOWARN_SIGNS
#include <boost/test/unit_test.hpp>
DYNDIST_NOWARN_POP

#include "dyndist/discrete_distribution.h"
#include "dyndist/utility.h"

#if __cplusplus >= 201103L

/* C++ 11. Use dyndist's RNG adapter for the standard library RNG */
#include "dyndist/rng_stdcxx.h"

#else // __cplusplus < 201103L

/* C++ 98. Use dyndist's RNG adapter for the boost RNG */
#include "dyndist/rng_boost.h"

#endif // __cplusplus <> 201103L

/**
 * \struct weight_flavour
 *
 * Dummy value used as "flavour" by flavoured_weight below.
 */
struct weight_flavour
{};

/*
 * Forward declarations
 */
template<typename ValueType>
struct flavoured_sample;

/**
 * \struct flavoured_weight
 *
 * Used to test discrete_distribution's support for "flavoured weighte", i.e.
 * for a WeightType for which only those instances are compatible which were
 * produced from one another by WeightType's operations.
 */
template<typename ValueType>
struct flavoured_weight
{
    typedef ValueType value_type;

    flavoured_weight()
        :value(0), flavour(NULL)
    {}

    flavoured_weight(const weight_flavour* _flavour)
        :value(0), flavour(_flavour)
    {}

    flavoured_weight(ValueType _value, const weight_flavour* _flavour)
        :value(_value), flavour(_flavour)
    {}
    
    flavoured_weight& operator=(value_type _value)
    { DYNDIST_ASSERT(_value == 0); value = _value; return *this; }

    flavoured_weight& operator+=(const flavoured_weight& o)
    { DYNDIST_ASSERT(flavour == o.flavour); value += o.value; return *this; }

    flavoured_weight& operator-=(const flavoured_weight& o)
    { DYNDIST_ASSERT(flavour == o.flavour); value -= o.value; return *this; }

    flavoured_weight operator+(flavoured_weight o) const
    { o += *this; return *this; }

    flavoured_weight operator-(flavoured_weight o) const
    { o -= *this; return *this; }

    bool operator< (const flavoured_weight& o) const
    { DYNDIST_ASSERT(flavour == o.flavour); return value <  o.value; }
    bool operator<=(const flavoured_weight& o) const
    { DYNDIST_ASSERT(flavour == o.flavour); return value <= o.value; }
    bool operator==(const flavoured_weight& o) const
    { DYNDIST_ASSERT(flavour == o.flavour); return value == o.value; }
    bool operator!=(const flavoured_weight& o) const
    { DYNDIST_ASSERT(flavour == o.flavour); return value != o.value; }
    bool operator>=(const flavoured_weight& o) const
    { DYNDIST_ASSERT(flavour == o.flavour); return value >= o.value; }
    bool operator> (const flavoured_weight& o) const
    { DYNDIST_ASSERT(flavour == o.flavour); return value >  o.value; }

    friend
    std::ptrdiff_t
    log2ceil(const flavoured_weight& weight)
    { return dyndist::log2ceil(weight.value); }

    friend
    flavoured_weight
    pow2(const std::ptrdiff_t exponent, const flavoured_weight& result_template)
    { return flavoured_weight(dyndist::pow2(exponent, result_template.value),
                              result_template.flavour); }

private:
    friend struct flavoured_sample<value_type>;
    
    template<typename,typename,typename>
    friend struct DYNDIST_NAMESPACE::uniform_distribution;
    
    value_type value;
    const weight_flavour* flavour;
};

/**
 * \struct flavoured_sample
 *
 * Result type used by `uniform_distribution<flavoured_weight<...>, Engine>`.
 *
 * Used to verify that `discrete_distribution` only requires the two operations
 * `result -= weight` and `result < weight`.
 */
template<typename ValueType>
struct flavoured_sample
{
    typedef ValueType value_type;
    
    bool operator< (const flavoured_weight<ValueType>& w) const
    { DYNDIST_ASSERT(flavour == w.flavour); return value < w.value; }
    
    flavoured_sample& operator-=(const flavoured_weight<ValueType>& w)
    {
        DYNDIST_ASSERT(flavour == w.flavour);
        value -= w.value;
        DYNDIST_ASSERT(value >= 0);
        return *this;
    }
    
private:
    template<typename,typename,typename>
    friend struct DYNDIST_NAMESPACE::uniform_distribution;
    
    flavoured_sample(ValueType _value, const weight_flavour* _flavour)
        :value(_value), flavour(_flavour)
    {}
    
    value_type value;
    const weight_flavour* flavour;
};

DYNDIST_NAMESPACE_BEGIN

/**
 * \struct uniform_distribution
 *
 * `uniform_distribution` implementation for `flavoured_weight<V>`. Simply
 * forwards to `V`'s `uniform_distribution` implementation.
 */
template<typename V, typename Engine>
struct uniform_distribution< flavoured_weight<V>, Engine >
{
    typedef flavoured_sample<V> result_type;
    
    static result_type sample(const flavoured_weight<V>& upper, Engine& engine)
    {
        const V v = uniform_distribution<V,Engine>::sample(upper.value, engine);
        return flavoured_sample<V>(v, upper.flavour);
    }
};

DYNDIST_NAMESPACE_END

template<typename WeightType>
struct vector_distribution_types
{
    typedef dyndist::discrete_distribution_pointer
    <std::size_t, WeightType, std::size_t>
    event_pointer;
    
    typedef std::vector<event_pointer> events_type;

    struct event_moved_callback
    {
        event_moved_callback(events_type& vector)
            :m_vector(vector)
        {}
        
        void operator()(event_pointer event, event_pointer event_new)
        {
            DYNDIST_ASSERT(event->data == event_new->data);
            DYNDIST_ASSERT(event->data < m_vector.size());
            DYNDIST_ASSERT(m_vector[event->data] == event);
            m_vector[event->data] = event_new;
        }
        
    private:
        events_type& m_vector;
    };
    
    typedef
    dyndist::discrete_distribution<std::size_t, WeightType, std::size_t, event_moved_callback>
    distribution_type;
};

template<typename WeightType>
class vector_distribution
    :public vector_distribution_types<WeightType>::distribution_type
{
    typedef typename vector_distribution_types<WeightType>::distribution_type super_type;
    
    typedef typename vector_distribution_types<WeightType>::events_type events_type;
    
public:
    typedef typename super_type::value_type value_type;
    
    typedef typename super_type::pointer pointer;
    
    typedef typename super_type::weight_type weight_type;
    
    vector_distribution()
        :super_type(typename vector_distribution_types<WeightType>::
                    event_moved_callback(m_events))
    {}
    
    void clear()
    { super_type::clear(); m_events.clear(); }
    
    void set(std::size_t index, const weight_type& weight);
    
    void erase(std::size_t index);
    
    template<typename Engine>
    std::size_t operator()(Engine& engine);
    
private:
    void verify_consistency();

    /* Hide */
    void insert();
    
    events_type m_events;
};

template<typename W>
void
vector_distribution<W>::set(std::size_t index, const weight_type& weight)
{
    if (index >= m_events.size())
        m_events.resize(index + 1);
    
    if (!m_events[index]) {
        const pointer p = super_type::insert(value_type(weight, index));
        DYNDIST_ASSERT(p->data == index);
        m_events[index] = p;
    }
    else {
        DYNDIST_ASSERT(m_events[index]->data == index);
        const pointer p = super_type::set(m_events[index], weight);
        DYNDIST_ASSERT(p->data == index);
        DYNDIST_ASSERT(m_events[index] == p);
    }

    verify_consistency();
}

template<typename W>
void
vector_distribution<W>::erase(std::size_t index)
{
    if (index >= m_events.size())
        return;
    
    if (m_events[index]) {
        super_type::erase(m_events[index]);
        m_events[index] = pointer();
    }
    
    while (!m_events.empty() && !m_events.back())
        m_events.pop_back();
    
    verify_consistency();
}

template<typename W>
template<typename Engine>
std::size_t
vector_distribution<W>::operator()(Engine& engine)
{
    const pointer s = super_type::operator()(engine);
    DYNDIST_ASSERT(m_events[s->data] == s);
    verify_consistency();
    
    return s->data;
}

template<typename W>
void
vector_distribution<W>::verify_consistency()
{
    super_type::verify_consistency();
    
    std::vector<bool> event_found = std::vector<bool>(m_events.size(), false);
    
    {
        const typename super_type::iterator end = super_type::end();
        for(typename super_type::iterator i = super_type::begin();
            i != end;
            ++i)
        {
            DYNDIST_ASSERT(i->data < m_events.size());
            DYNDIST_ASSERT(m_events[i->data] == &*i);
            DYNDIST_ASSERT(*m_events[i->data] == *i);
            DYNDIST_ASSERT(!event_found[i->data]);
            event_found[i->data] = true;
        }
    }
    
    {
        const typename events_type::const_iterator event_end = m_events.end();
        const typename std::vector<bool>::const_iterator found_end = event_found.end();
        typename events_type::const_iterator event_i;
        typename std::vector<bool>::const_iterator found_i;
        for(event_i = m_events.begin(), found_i = event_found.begin();
            (event_i != event_end) && (found_i != found_end);
            ++event_i, ++found_i)
        {
            DYNDIST_ASSERT(*found_i == (bool)*event_i);
        }
    }
}

BOOST_AUTO_TEST_SUITE(test_discrete_distribution)

namespace {
    template<typename D>
    void
    test_iteration(D& dist)
    {
        typedef typename D::iterator iterator;
        typedef typename D::pointer pointer;
        typedef typename D::size_t size_t;
        size_t n = 0;
        for(iterator e = dist.begin(); e != dist.end(); ++e, ++n) {
            BOOST_REQUIRE(e == dist.begin() + n);
            BOOST_REQUIRE(e == dist.end() - (dist.size() - n));
            BOOST_REQUIRE(e - dist.begin() == n);
            BOOST_REQUIRE(dist.begin() - e == -(std::ptrdiff_t)n);
            BOOST_REQUIRE(dist.end() - e == dist.size() - n);
            BOOST_REQUIRE(e - dist.end() == -(std::ptrdiff_t)(dist.size() - n));
            
            iterator r = e;
            for(size_t j=n; j > 0; --j, --r) {
                BOOST_REQUIRE(r + (n - j) == e);
                BOOST_REQUIRE(r - e == -(std::ptrdiff_t)(n - j));
                BOOST_REQUIRE(r <= e);
                BOOST_REQUIRE((j == n) || (r < e));
            }
            BOOST_REQUIRE(r == dist.begin());
            
            for(size_t j=0; j < n; ++j, ++r) {
                BOOST_REQUIRE(e - (n - j) == r);
                BOOST_REQUIRE(e - r == n - j);
                BOOST_REQUIRE(!(r >= e));
                BOOST_REQUIRE(!(r > e));
            }
            BOOST_CHECK(r == e);
            
            for(size_t j=n; j < dist.size(); ++j, ++r) {
                BOOST_REQUIRE(r + -(std::ptrdiff_t)(j - n) == e);
                BOOST_REQUIRE(r - e == j - n);
                BOOST_REQUIRE(e <= r);
                BOOST_REQUIRE((j == n) || (e < r));
            }
            BOOST_REQUIRE(r == dist.end());
            
            BOOST_REQUIRE(&*e == &dist.begin()[n]);
            const pointer p = e;
            BOOST_REQUIRE(&*p == &*e);
        }
    }
    
    void test_empty_distribution(bool use_meta)
    {
        typedef flavoured_weight<unsigned int> weight_type;
        weight_flavour flavour;
        vector_distribution<weight_type> dist;
        dist.force_meta(use_meta ? 1 : -1);
        
        BOOST_CHECK(dist.weight() == weight_type());
        BOOST_CHECK(dist.begin() == dist.end());
        
        {
            bool did_throw_out_of_range = false;
            try {
                mt19937 rng;
                dist(rng);
            }
            catch (const std::out_of_range&) {
                did_throw_out_of_range = true;
            }
            BOOST_REQUIRE(did_throw_out_of_range);
        }
        
        dist.set(0, weight_type(0, &flavour));
        BOOST_CHECK(dist.weight() == weight_type(&flavour));
        BOOST_CHECK(dist.size() == 1);
        BOOST_CHECK(dist.size_nonzero() == 0);
        BOOST_CHECK(++dist.begin() == dist.end());
        test_iteration(dist);
        
        {
            bool did_throw_out_of_range = false;
            try {
                mt19937 rng;
                dist(rng);
            }
            catch (const std::out_of_range&) {
                did_throw_out_of_range = true;
            }
            BOOST_REQUIRE(did_throw_out_of_range);
        }
        
        dist.erase(0);
        BOOST_CHECK(dist.weight() == weight_type(&flavour));
        BOOST_CHECK(dist.size() == 0);
        BOOST_CHECK(dist.size_nonzero() == 0);
        BOOST_CHECK(dist.begin() == dist.end());
    }
}

BOOST_AUTO_TEST_CASE(empty_distribution_nometa)
{
    test_empty_distribution(false);
}

BOOST_AUTO_TEST_CASE(empty_distribution_meta)
{
    test_empty_distribution(true);
}

namespace {
    void test_set(bool use_meta)
    {
        mt19937 rng;
        vector_distribution<unsigned int> dist;
        dist.force_meta(use_meta ? 1 : -1);

        test_iteration(dist);
        
        dist.set(0, 1); BOOST_CHECK_EQUAL(dist.size(), 1); BOOST_CHECK_EQUAL(dist.weight(), 1);
        test_iteration(dist);
        dist(rng);
        
        dist.set(1, 1); BOOST_CHECK_EQUAL(dist.size(), 2); BOOST_CHECK_EQUAL(dist.weight(), 2);
        test_iteration(dist);
        dist(rng);
        
        dist.set(2, 1); BOOST_CHECK_EQUAL(dist.size(), 3); BOOST_CHECK_EQUAL(dist.weight(), 3);
        test_iteration(dist);
        dist(rng);
        
        dist.set(0, 2); BOOST_CHECK_EQUAL(dist.weight(), 4);
        test_iteration(dist);
        dist(rng);
        
        dist.set(1, 1); BOOST_CHECK_EQUAL(dist.weight(), 4);
        test_iteration(dist);
        dist(rng);
        
        dist.set(2, 0); BOOST_CHECK_EQUAL(dist.weight(), 3);
        test_iteration(dist);
        dist(rng);
        
        dist.set(0, 0); BOOST_CHECK_EQUAL(dist.weight(), 1);
        test_iteration(dist);
        dist(rng);
        
        dist.set(1, 0); BOOST_CHECK_EQUAL(dist.weight(), 0);
        test_iteration(dist);
        BOOST_CHECK_EQUAL(dist.size_nonzero(), 0);
        
        dist.set(2, 0); BOOST_CHECK_EQUAL(dist.weight(), 0);
        test_iteration(dist);
        BOOST_CHECK_EQUAL(dist.size_nonzero(), 0);
        
        dist.erase(1); BOOST_CHECK_EQUAL(dist.size(), 2);
        test_iteration(dist);
        
        dist.erase(2); BOOST_CHECK_EQUAL(dist.size(), 1);
        test_iteration(dist);
        
        dist.erase(0); BOOST_CHECK_EQUAL(dist.size(), 0);
        test_iteration(dist);
        
        dist.set(0, 1); BOOST_CHECK_EQUAL(dist.size(), 1); BOOST_CHECK_EQUAL(dist.weight(), 1);
        test_iteration(dist);
        dist(rng);
        
        dist.set(1, 1); BOOST_CHECK_EQUAL(dist.size(), 2); BOOST_CHECK_EQUAL(dist.weight(), 2);
        test_iteration(dist);
        dist(rng);
        
        dist.set(2, 1); BOOST_CHECK_EQUAL(dist.size(), 3); BOOST_CHECK_EQUAL(dist.weight(), 3);
        test_iteration(dist);
        dist(rng);
    }
}

BOOST_AUTO_TEST_CASE(set_nometa)
{
    test_set(false);
}

BOOST_AUTO_TEST_CASE(set_meta)
{
    test_set(true);
}

namespace {
    void test_draw_purify(bool use_meta)
    {
        using namespace dyndist;
        
        mt19937 rng;
        vector_distribution<uint64_t> dist;
        dist.force_meta(use_meta ? 1 : -1);
        
        dist.set(0, 256);
        dist.set(1, 2);
        dist.set(2, 2);
        dist.set(3, 16);
        dist.set(1, 1);
        dist.set(2, 1);
        dist.set(0, 0);
        
        /* Run until we draw 1 or 2 */
        while (true) {
            std::size_t s = dist(rng);
            BOOST_CHECK((s > 0) && (s <= 3));
            if (s != 3)
                break;
            if (use_meta) {
                BOOST_CHECK_EQUAL(dist.stats_draws_level_scan(), 0);
                BOOST_CHECK_EQUAL(dist.stats_draws_level_meta(), dist.stats_draws());
            }
            else {
                BOOST_CHECK_EQUAL(dist.stats_draws_level_scan(), dist.stats_draws());
                BOOST_CHECK_EQUAL(dist.stats_draws_level_meta(), 0);
            }
            BOOST_CHECK_EQUAL(dist.stats_draws_acceptor(), 0);
            BOOST_CHECK_EQUAL(dist.stats_draws_single(), dist.stats_draws());
            BOOST_CHECK_EQUAL(dist.stats_draws_scan(), 0);
            BOOST_CHECK_EQUAL(dist.stats_draws_purify(), 0);
        }
        
        /* Should have purified the level containing 1,2 now */
        if (use_meta) {
            BOOST_CHECK_EQUAL(dist.stats_draws_level_scan(), 0);
            BOOST_CHECK_EQUAL(dist.stats_draws_level_meta(), dist.stats_draws());
        }
        else {
            BOOST_CHECK_EQUAL(dist.stats_draws_level_scan(), dist.stats_draws());
            BOOST_CHECK_EQUAL(dist.stats_draws_level_meta(), 0);
        }
        BOOST_CHECK_EQUAL(dist.stats_draws_acceptor(), 0);
        BOOST_CHECK_EQUAL(dist.stats_draws_single(), dist.stats_draws() - 1);
        BOOST_CHECK_EQUAL(dist.stats_draws_scan(), 0);
        BOOST_CHECK_EQUAL(dist.stats_draws_purify(), 1);
        
        /* Further draws shouldn't purify again */
        const std::size_t N = 100;
        for(std::size_t i = 0; i < N; ++i) {
            const std::size_t s = dist(rng);
            BOOST_CHECK((s > 0) && (s <= 3));
        }
        if (use_meta) {
            BOOST_CHECK_EQUAL(dist.stats_draws_level_scan(), 0);
            BOOST_CHECK_EQUAL(dist.stats_draws_level_meta(), dist.stats_draws());
        }
        else {
            BOOST_CHECK_EQUAL(dist.stats_draws_level_scan(), dist.stats_draws());
            BOOST_CHECK_EQUAL(dist.stats_draws_level_meta(), 0);
        }
        BOOST_CHECK_EQUAL(dist.stats_draws_acceptor() +
                          dist.stats_draws_single() +
                          dist.stats_draws_scan(),
                          dist.stats_draws() - 1);
        BOOST_CHECK_EQUAL(dist.stats_draws_purify(), 1);
    }
}

BOOST_AUTO_TEST_CASE(draw_purify_nometa)
{
    test_draw_purify(false);
}

BOOST_AUTO_TEST_CASE(draw_purify_meta)
{
    test_draw_purify(true);
}

namespace {
    void test_draw_impure(bool use_meta)
    {
        using namespace dyndist;
        
        mt19937 rng;
        vector_distribution<std::size_t> dist;
        dist.force_meta(use_meta ? 1 : -1);
        
        dist.set(0, 256);
        dist.set(1, 2);
        dist.set(2, 2);
        dist.set(1, 1);
        dist.set(2, 1);
        while (dist(rng) == 0);
        if (use_meta) {
            BOOST_CHECK_EQUAL(dist.stats_draws_level_scan(), 0);
            BOOST_CHECK_EQUAL(dist.stats_draws_level_meta(), dist.stats_draws());
        }
        else {
            BOOST_CHECK_EQUAL(dist.stats_draws_level_scan(), dist.stats_draws());
            BOOST_CHECK_EQUAL(dist.stats_draws_level_meta(), 0);
        }
        BOOST_CHECK_EQUAL(dist.stats_draws_single(), dist.stats_draws() - 1);
        BOOST_CHECK_EQUAL(dist.stats_draws_acceptor(), 0);
        BOOST_CHECK_EQUAL(dist.stats_draws_scan(), 1);
        BOOST_CHECK_EQUAL(dist.stats_draws_purify(), 0);
    }
}

BOOST_AUTO_TEST_CASE(draw_impure_nometa)
{
    test_draw_impure(false);
}

BOOST_AUTO_TEST_CASE(draw_impure_meta)
{
    test_draw_impure(true);
}

BOOST_AUTO_TEST_CASE(set_random_uniform_nometa)
{
    using namespace dyndist;
    
    vector_distribution<uint64_t> dist;
    dist.force_meta(-1);
    
    const std::size_t N = 1000;
    mt19937 rng;
    uniform_int_distribution<uint64_t> dist_key(0,N-1);
    const uint64_t max_weight = std::numeric_limits<uint64_t>::max() / N;
    uniform_int_distribution<uint64_t> dist_weight(0, max_weight);
    
    for(std::size_t i=0; i < N; ++i) {
        dist.set(dist_key(rng), dist_weight(rng));
        const uint64_t dr = dist(rng);
        BOOST_REQUIRE(0 <= dr);
        BOOST_REQUIRE(dr < N);
    }
    
    test_iteration(dist);
}

BOOST_AUTO_TEST_CASE(set_random_uniform_meta)
{
    using namespace dyndist;
    
    vector_distribution<uint64_t> dist;
    dist.force_meta(1);
    
    const std::size_t N = 1000;
    mt19937 rng;
    uniform_int_distribution<uint64_t> dist_key(0,N-1);
    const uint64_t max_weight = std::numeric_limits<uint64_t>::max() / N;
    uniform_int_distribution<uint64_t> dist_weight(0, max_weight);
    
    for(std::size_t i=0; i < N; ++i) {
        dist.set(dist_key(rng), dist_weight(rng));
        const uint64_t dr = dist(rng);
        BOOST_REQUIRE(0 <= dr);
        BOOST_REQUIRE(dr < N);
    }
    
    test_iteration(dist);
}

BOOST_AUTO_TEST_CASE(set_random_loguniform_nometa)
{
    using namespace dyndist;
    
    vector_distribution<uint64_t> dist;
    dist.force_meta(-1);
    
    const std::size_t N = 1000;
    mt19937 rng;
    uniform_int_distribution<uint64_t> dist_key(0,N-1);
    const uint64_t max_weight = std::numeric_limits<uint64_t>::max() / N;
    uniform_real_distribution<double> dist_logweight(0.0, std::log(max_weight));
    
    for(std::size_t i=0; i < N; ++i) {
        dist.set(dist_key(rng), (uint64_t)std::exp(dist_logweight(rng)));
        const uint64_t dr = dist(rng);
        BOOST_REQUIRE(0 <= dr);
        BOOST_REQUIRE(dr < N);
    }
    
    test_iteration(dist);
}

BOOST_AUTO_TEST_CASE(set_random_loguniform_meta)
{
    using namespace dyndist;
    
    vector_distribution<uint64_t> dist;
    dist.force_meta(1);
    
    const std::size_t N = 1000;
    mt19937 rng;
    uniform_int_distribution<uint64_t> dist_key(0,N-1);
    const uint64_t max_weight = std::numeric_limits<uint64_t>::max() / N;
    uniform_real_distribution<double> dist_logweight(0.0, std::log(max_weight));
    
    for(std::size_t i=0; i < N; ++i) {
        dist.set(dist_key(rng), (uint64_t)std::exp(dist_logweight(rng)));
        const uint64_t dr = dist(rng);
        BOOST_REQUIRE(0 <= dr);
        BOOST_REQUIRE(dr < N);
    }
    
    test_iteration(dist);
}

BOOST_AUTO_TEST_SUITE_END()
