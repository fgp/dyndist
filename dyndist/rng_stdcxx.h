/*
 * dyndist/rng_stdcxx.h
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

#ifndef dyndist_rng_boost_h
#define dyndist_rng_boost_h

#include <type_traits>
#include <random>

#include "dyndist/global.h"
#include "dyndist/utility.h"

DYNDIST_NAMESPACE_BEGIN

/**
 * uniform_distribution<floating_point_type, Engine>
 *
 * Uses boost's uniform_real_distribution to draw uniformy distributed
 * floating point values. uniform_real_distribution<T> must be valid, and
 * Engine must be one of boost's random number generators.
 */
template<typename T, typename Engine>
struct uniform_distribution
<T, Engine, typename std::enable_if<std::is_floating_point<T>::value>::type>
{
    typedef T result_type;

    DYNDIST_INLINE
    static result_type sample(const T upper, Engine& engine)
    {
        typedef std::uniform_real_distribution<T> dist_t;
        return dist_t((T)0, upper)(engine);
    }
};

/**
 * uniform_distribution<integral_type, Engine>
 *
 * Uses boost's uniform_int_distribution to draw uniformy distributed
 * integral values. uniform_int_distribution<T> must be valid, and
 * Engine must be one of boost's random number generators.
 */
template<typename T, typename Engine>
struct uniform_distribution
<T, Engine, typename std::enable_if<std::is_integral<T>::value>::type>
{
    typedef T result_type;

    DYNDIST_INLINE
    static result_type sample(const T upper, Engine& engine)
    {
        typedef std::uniform_int_distribution<T> dist_t;
        DYNDIST_ASSERT(upper >= T(1));
        return dist_t((T)0, upper - T(1))(engine);
    }
};

DYNDIST_NAMESPACE_END

#endif
