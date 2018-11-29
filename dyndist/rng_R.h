/*
 * dyndist/rng_R.h
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

#ifndef dyndist_rng_R_h
#define dyndist_rng_R_h

#include "dyndist/global.h"

#include <cmath>
#include <utility>
#include <limits>

#include <Rcpp.h>

#include "dyndist/utility.h"

DYNDIST_NAMESPACE_BEGIN

/**
 * rng_R
 *
 * Type representing the R random number generator. To be used as the Engine type
 * with uniform_distribution and discrete_distribution.
 */
struct rng_R {
    rng_R() {}
    
private:
    /* Shouldn't be copied */
    rng_R(const rng_R&);
    
    Rcpp::RNGScope rng_scope;
};

/**
 * uniform_distribution<double, rng_R>
 *
 * Specialization of dyndist::uniform_distribution which uses the R random number
 * generator to draw uniformly distributed double values.
 */
template<>
struct uniform_distribution<double, rng_R, void>
{
    typedef double weight_type;
    typedef rng_R engine_type;
    typedef weight_type bound_type;
    typedef weight_type result_type;
    
    DYNDIST_INLINE
    static result_type sample(const double upper, rng_R&)
    {
        DYNDIST_ASSERT(upper >= 0.0);
        return (unif_rand() * upper);
    }
};

/**
 * uniform_distribution<double, rng_R>
 *
 * Specialization of dyndist::uniform_distribution which uses the R random number
 * generator to draw uniformly integer values.
 */
template<>
struct uniform_distribution<int, rng_R, void>
{
    typedef int weight_type;
    typedef rng_R engine_type;
    typedef weight_type bound_type;
    typedef weight_type result_type;
    
    DYNDIST_INLINE
    static result_type sample(const int upper, rng_R&)
    {
        DYNDIST_ASSERT(upper >= 0);
        return (int)std::floor(unif_rand() * double(upper));
    }
};

/**
 * uniform_distribution<pair<int, pow2_expression<int>, rng_R>
 *
 * Specialization of dyndist::uniform_distribution which uses the R random number
 * generator to draw pairs of uniformly distributed integer values from
 * [0, b) x [0, 2^e). This is an optimization for discrete_distribution, which
 * requires such pairs fairly often, and it thus pays off to generate them using
 * only one call to the underlying RNG if possible. See the explanations in
 * discrete_distribution.h for details.
 */

template<>
struct uniform_distribution<std::pair< int, pow2_expression<int> >, rng_R, void>
{
    typedef rng_R engine_type;
    typedef std::pair< int, pow2_expression<int> > bound_type;
    typedef std::pair< int, int > result_type;

    static result_type sample(const bound_type& upper, rng_R& engine);
};

DYNDIST_INLINE
uniform_distribution<std::pair< int, pow2_expression<int> >, rng_R, void>
::result_type
uniform_distribution<std::pair< int, pow2_expression<int> >, rng_R, void>
::sample
(const bound_type& upper, rng_R& engine)
{
    typedef uniform_distribution<int, rng_R> unif_t;
    
    DYNDIST_ASSERT(0 < upper.first);
    DYNDIST_ASSERT(0 < upper.second);
    const std::ptrdiff_t e = upper.second.exponent;
    DYNDIST_ASSERT(e >= 0);
    DYNDIST_ASSERT((int)upper.second == pow2(e, (int)0));
    
    const unsigned int u = ((unsigned int)upper.first) << e;
    const unsigned int u_max = (unsigned int)std::numeric_limits<int>::max();
    const bool u_valid = ((u <= u_max) && ((((int)u) >> e) == upper.first ));
    if (u_valid) {
        /* Draw combined, u = u1 * u2 <= std::numeric_limits<int>::max() */
        DYNDIST_ASSERT(u / (int)upper.second == upper.first);
        DYNDIST_ASSERT(u % (int)upper.second == 0);
        
        /* Draw from [0, u1 * u2) */
        const int r = unif_t::sample((int)u, engine);
        DYNDIST_ASSERT((r >= 0) && (r < u));
        
        /* Return (r / u1, r % u1) */
        const int r1 = r >> e;
        DYNDIST_ASSERT(r1 < upper.first);
        const int r2 = r & ((((int)1) << e) - 1);
        DYNDIST_ASSERT(r2 < upper.second);
        return result_type(r1, r2);
    }
    else {
        /* Draw separately */
        const int r1 = unif_t::sample(upper.first, engine);
        const int r2 = unif_t::sample(upper.second, engine);
        return result_type(r1, r2);
    }
}

DYNDIST_NAMESPACE_END

#endif
