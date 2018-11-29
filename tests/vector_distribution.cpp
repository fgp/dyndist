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

BOOST_AUTO_TEST_SUITE(test_vector_distribution)

BOOST_AUTO_TEST_CASE(chisquared_test)
{
    using namespace dyndist;

    const double alpha = 1e-6;
    const std::size_t M = 10;
    const std::size_t N = 100;
    const std::size_t S = 1000000;
    mt19937 rng;
    uniform_int_distribution<uint64_t> dist_key(0,N-1);
    uniform_int_distribution<uint64_t> dist_weight(1, 1000);

    vector_distribution<uint64_t> p(N, 0);
    std::vector<uint64_t> obs;

    for(std::size_t m=0; m < M; ++m) {
        /* Re-assign weights */
        for (std::size_t i = 0; i < N; ++i) {
            p[dist_key(rng)] = dist_weight(rng);
        }

        /* Samples B times */
        obs = std::vector<uint64_t>(N, 0);
        for (std::size_t s = 0; s < S; ++s)
            obs[p(rng)] += 1;

        /* Check that the empirical and expected distributions are similar */
        double x = 0;
        int df = -1;
        for (std::size_t i = 0; i < N; ++i) {
            if (p[i] == 0)
                continue;

            const double m = (double) p[i] * S / p.weight();
            x += std::pow((double)obs[i] - m, 2) / m;
            ++df;
        }

        /* Check that the p-value is reasonable */
        boost::math::chi_squared chisq_dist(df);
        const double critical = quantile(complement(chisq_dist, alpha));
        BOOST_REQUIRE_LT(x, critical);
    }
}

BOOST_AUTO_TEST_SUITE_END()
