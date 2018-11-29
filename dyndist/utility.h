/*
 * dyndist/utility.h
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
#ifndef dyndist_utility_h
#define dyndist_utility_h

#include "dyndist/global.h"

#include <cmath>
#include <utility>

DYNDIST_NAMESPACE_BEGIN

/**
 * \fn sum_out_of_range
 *
 * Returns true if the sum of two values isn't representable by type of the two values
 */
template <typename IntegralType>
bool sum_out_of_range(const IntegralType a, const IntegralType b) {
    /* See: https://www.securecoding.cert.org/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow */
    if ((b > 0) && (a > (std::numeric_limits<IntegralType>::max() - b)))
        return true;
    if ((b < 0) && (a < (std::numeric_limits<IntegralType>::min() - b)))
        return true;
    return false;
}

/**
 * \fn diff_out_of_range
 *
 * Returns true if the difference of two values isn't representable by type of the two values
 */
template <typename IntegralType>
bool diff_out_of_range(const IntegralType a, const IntegralType b) {
    /* See: https://www.securecoding.cert.org/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow */
    if ((b > 0) && (a < (std::numeric_limits<IntegralType>::min() + b)))
        return true;
    if ((b < 0) && (a > (std::numeric_limits<IntegralType>::max() + b)))
        return true;
    return false;
}

/**
 * \fn log2ceil
 *
 * Generic log2ceil implementation, which computes ceil(log2(value)) for all
 * types which can be converted to double.
 */
template<typename WeightType>
DYNDIST_INLINE
std::ptrdiff_t
log2ceil(const WeightType& value)
{
    /* XXX: Rounding behaviour? */
    int exponent = 0;
    const double mantissa = std::frexp((double)value, &exponent);
    return exponent - ((mantissa == 0.5L) ? 1 : 0);
}

/**
 * \fn pow2
 *
 * Generic pow2 implementation, which computes V(2) ^ exponent for all
 * integral types which provide a left-shift operator.
 */
template<typename WeightType>
DYNDIST_INLINE
WeightType
pow2(const std::ptrdiff_t exponent, const WeightType&)
{
    DYNDIST_ASSERT(exponent < std::numeric_limits<WeightType>::digits);
    return (exponent >= 0) ? (WeightType(1) << exponent) : WeightType(0);
}

/**
 * \fn pow2_expression
 *
 * Represents an expression of the form 2^exponent. Can be used to optimize
 * cases where the argument to a function is always an integral power of two,
 * and where evaluation of the function is much cheaper for such arguments.
 *
 * Currently only used in conjunction with uniform_distribution, and then only
 * to draw pairs of random numbers, where one of the upper bounds is a power of
 * 2. See the uniform_distribution specialization below.
 */
template<typename WeightType>
struct pow2_expression
{
    typedef WeightType weight_type;
    
    pow2_expression(const std::ptrdiff_t _exponent,
                    const weight_type& _result)
        :exponent(_exponent), result(_result)
    {
        DYNDIST_ASSERT(result == pow2(exponent, result));
    }
    
    DYNDIST_INLINE
    operator const weight_type& () const
    { return result; }
    
    const std::ptrdiff_t exponent;
    const weight_type result;
};

/**
 * \struct uniform_distribution
 *
 * Generic uniform_distribution dummy. Must be specialized for specific
 * combinations of WeightType and Engine. The sample() method of this generic
 * dummy implementation isn't defined, only declared, meaning that missing
 * specializations will cause linker errors.
 *
 * Note that `uniform_distribution` types are never instantiated -- the
 * `sample()` member is static!
 */
template<typename WeightType, typename Engine, typename Enable = void>
struct uniform_distribution
{
    typedef WeightType weight_type;
    typedef Engine engine_type;
    typedef weight_type bound_type;
    typedef weight_type result_type;
    
    static result_type sample(const bound_type& upper, engine_type& engine);
    
private:
    uniform_distribution();
};

/**
 * \struct uniform_distribution< pow2_expression<WeightType> >
 *
 * Generic implementation of uniform_distribution for bounds which are
 * integral powers of two. Simply forwards to the specialization for arbitrary
 * bounds of type WeightType. See discrete_distribution for why and when it
 * makes sense to explicitly specialize this case.
 */
template<typename WeightType, typename Engine>
struct uniform_distribution< pow2_expression<WeightType>, Engine>
{
private:
    typedef uniform_distribution<WeightType, Engine> unif_type;
    
    uniform_distribution();
    
public:
    typedef WeightType weight_type;
    typedef Engine engine_type;
    typedef pow2_expression<weight_type> bound_type;
    typedef typename unif_type::result_type result_type;
    
    DYNDIST_INLINE
    static result_type sample(const bound_type& upper, Engine& engine)
    {
        return unif_type::sample(upper, engine);
    }
};

/**
 * \struct uniform_distribution< std::pair<WeightType1, WeightType2> >
 *
 * Generic implementation of uniform_distribution for pairs of random numbers.
 * Simply forwards to the specializations for the individual weigh types,
 * and combines their results into a pair. See discrete_distribution for why
 * and when it makes sense to explicitly specialize this case.
 */
template<typename WeightType1, typename WeightType2, typename Engine>
struct uniform_distribution< std::pair<WeightType1, WeightType2>, Engine >
{
private:
    typedef uniform_distribution<WeightType1, Engine> unif1_type;
    typedef uniform_distribution<WeightType2, Engine> unif2_type;
    typedef typename unif1_type::result_type result1_type;
    typedef typename unif2_type::result_type result2_type;
    
    uniform_distribution();
    
public:
    typedef Engine engine_type;
    typedef std::pair<WeightType1, WeightType2> bound_type;
    typedef std::pair<result1_type, result2_type> result_type;
    
    DYNDIST_INLINE
    static result_type sample(const bound_type& upper, Engine& engine)
    {
        return result_type(unif1_type::sample(upper.first, engine),
                           unif2_type::sample(upper.second, engine));
    }
};

DYNDIST_NAMESPACE_END

#endif
