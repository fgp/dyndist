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

#include "dyndist/global.h"

#include <cmath>
#include <limits>
#include <iostream>

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
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/test/unit_test.hpp>
DYNDIST_NOWARN_POP

#include "dyndist/vector_distribution.h"
#include "dyndist/utility.h"

#if __cplusplus >= 201103L

/* C++ 11. Use dyndist's RNG adapter for the standard library RNG */
#include "dyndist/rng_stdcxx.h"

#else // __cplusplus < 201103L

/* C++ 98. Use dyndist's RNG adapter for the boost RNG */
#include "dyndist/rng_boost.h"

#endif // __cplusplus <> 201103L

/* Instantiate templates to get correct code coverage results */
DYNDIST_NAMESPACE_BEGIN
template class vector_distribution<std::size_t>;
DYNDIST_NAMESPACE_END

BOOST_AUTO_TEST_SUITE(test_vector_distribution)

BOOST_AUTO_TEST_CASE(iterator_test)
{
    using namespace dyndist;

    typedef vector_distribution <uint64_t> dist_t;
    dist_t d;

    BOOST_REQUIRE(d.begin() == d.end());
    BOOST_REQUIRE(d.cbegin() == d.cend());
    //BOOST_REQUIRE(d.begin() == d.cend()); // currently unsupported
    BOOST_REQUIRE(d.cbegin() == d.end());

    const std::size_t N = 1000;
    for(std::size_t n = 0; n < N; ++n) {
        d.push_back(n);
        BOOST_REQUIRE_EQUAL(d.end() - d.begin(), n+1);
        BOOST_REQUIRE_EQUAL(d.cend() - d.cbegin(), n+1);
        //BOOST_REQUIRE_EQUAL(d.end() - d.cbegin(), n+1); // currently unsupported
        BOOST_REQUIRE_EQUAL(d.cend() - d.begin(), n+1);

        BOOST_REQUIRE(++(d.begin() + n) == d.end());
        BOOST_REQUIRE(++(d.cbegin() + n) == d.cend());

        BOOST_REQUIRE(--(d.end() - n) == d.begin());
        BOOST_REQUIRE(--(d.cend() - n) == d.cbegin());

        dist_t::iterator i1 = d.begin(); i1 += n; i1++; BOOST_REQUIRE(i1 == d.end());
        dist_t::const_iterator i2 = d.cbegin(); i2 += n; i2++; BOOST_REQUIRE(i2 == d.cend());

        BOOST_REQUIRE(d.begin() < d.end());
        BOOST_REQUIRE(d.begin() <= d.end());
        BOOST_REQUIRE(d.end() >= d.begin());
        BOOST_REQUIRE(d.end() > d.begin());
        BOOST_REQUIRE(d.cbegin() < d.cend());
        BOOST_REQUIRE(d.cbegin() <= d.cend());
        BOOST_REQUIRE(d.cend() >= d.cbegin());
        BOOST_REQUIRE(d.cend() > d.cbegin());

        BOOST_REQUIRE(&d[n] + 1 == d.end());
    }
}

BOOST_AUTO_TEST_CASE(chisquared_test)
{
    using namespace dyndist;

    const double alpha = 1e-3;
    const std::size_t R = 5;
    const std::size_t M = 10;
    const std::size_t N = 100;
    const std::size_t S = 100000;
    mt19937 rng;
    uniform_int_distribution<uint64_t> dist_key(0,N-1);
    uniform_int_distribution<uint64_t> dist_weight(1, 1000);

    for(std::size_t r=0; r < R; ++r) {
        typedef vector_distribution<uint64_t> dist_t;
        dist_t p(N, 0);
        std::vector<uint64_t> obs;

        double pval_max = 0.0;
        for(std::size_t m=0; m < M; ++m) {
            /* Re-assign weights */
            for (dist_t::iterator p_i = p.begin(); p_i != p.end(); ++p_i) {
                *p_i = dist_weight(rng);
            }

            /* Samples B times */
            obs = std::vector<uint64_t>(N, 0);
            for (std::size_t s = 0; s < S; ++s)
                obs[p(rng)] += 1;

            /* Check that the empirical and expected distributions are similar */
            double x = 0;
            int df = -1;
            for (dist_t::const_iterator p_i = p.cbegin(), p_end=p.cend(); p_i != p_end; ++p_i) {
                if (*p_i == 0)
                    continue;

                const double m = (double) *p_i * S / p.weight();
                x += std::pow((double)obs[p_i - p.cbegin()] - m, 2) / m;
                ++df;
            }

            /* Check that the p-value is reasonable */
            boost::math::chi_squared chisq_dist(df);
            const double pval = cdf(complement(chisq_dist, x));
            pval_max = std::max(pval, pval_max);
        }

        /* We should have observed at least one large p-value */
        BOOST_REQUIRE_GT(pval_max, 0.5);
    }
}

BOOST_AUTO_TEST_SUITE_END()
