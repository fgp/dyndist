/*
 * dyndist/discrete_distribution.h
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

/**
 * \class discrete_distribution
 *
 * Synopsis:
 * =========
 * discrete_distribution models an arbitrary discrete probability distribution
 * over an arbitrary finite set of values. Each value represents an event, and
 * has an associated weight as well as optional additional data. Weights must be
 * non-negative, but don't necessarily have to sum up to 1 -- i.e., the actual
 * probability of a value is the quotient of its weight and the sum of all the
 * value's weights.
 *
 * Interface:
 * ==========
 *
 * ### Type instantiation:
 * If the lifespan of references to values contained within a distribution never
 * extend past the invocation of non-const members (except `begin()` and `end()`,
 * that is), specifying an `EventMovedFunctor` will generally not be necessary,
 * and a concrete distribution type can simply be instatiated with
 *
 *     discrete_distribution<SizeType, WeightType, DataType>
 *
 * If the lifespan of references extends past calls of non-const members, users
 * of `discrete_distribution` will need to specify an `EventMovedFunctor` to be
 * notified if such invalidations occur (see the section on event pointers for
 * details). To be able to define a suitable ``EventMovedFunctor` type before
 * a concrete `discrete_distribution` type is instantiated, the distribution's
 * pointer type is a separate template, and can be instantiated with
 *
 *     typedef discrete_distribution_pointer<SizeType, WeightType, DataType>
 *             pointer_t;
 *
 * Using that type, a suitable `EventMovedFunctor` can be declared, usually
 * either as a function pointer, i.e.
 *
 *     typedef void (*ev_moved_t)(pointer_t ptr_old, pointer_t ptr_new);
 *
 * or as a type with an `operator()`, i.e.
 *
 *     struct ev_moved_t {
 *         void operator()(pointer_t ptr_old, pointer_t ptr_new);
 *     }
 *
 * A concrete distribution type is then instantiated with
 *
 *     typedef discrete_distribution<SizeType, WeightType, DataType, ev_moved_t>
 *             dist_t;
 *
 * ### Instantiation:
 *
 *     dist_t dist(ev_moved_t(...));
 *
 * or
 *
 *     dist_t dist();
 *
 * if `ev_moved_t` is default constructible.
 *
 * ### Insert event:
 *
 *     event_pointer e = dist.insert(dist_t::value_type(weight, data));
 *
 * ### Erase event:
 *
 *     dist.erase(e);
 *
 * ### Update weight:
 *
 *     dist.set(e, new_weight);
 *
 * ### Query number of events:
 *
 *     dist.size()
 *
 * ### Query number of events with non-zero weight:
 *
 *     dist.size_nonzero()
 *
 * ### Query sum of event weights:
 *
 *     dist.weight()
 *
 * ### Query existance of events with non-zero weight:
 *
 *     dist.empty()
 *
 * ### Drawing random events (throws if dist.empty()):
 *
 *     Engine engine;
 *     event_pointer random = dist(engine);
 *
 * ### Drawing random events, passing a random number (throws if dist.empty()):
 *
 *     typedef uniform_distribution<WeightType, Engine> rng_t;
 *     typedef rng_t::result_type random_t;
 *     Engine engine;
 *     random_t random = rng_t::sample(dist.weight(), engine);
 *     event_pointer random = dist(engine, random);
 *
 * Note that the distribution may still need to produce additional random
 * numbers!
 *
 * ### Iteration:
 *
 *     typedef dist_t::iterator iterator_t;
 *     const iterator_t end = dist.end();
 *     for(iterator_t i = dist.begin(); i != end; ++i) {
 *       const WeightType& i_weight = i->weight;
 *       DataType& i_data = i->data;
 *       pointer_t e = i; // but the reverse doesn't work
 *       value_t& v = *i ; // ... = *e would yield the same result
 *       ...
 *     }
 *
 * Iterators are allow random access, but the complexity of `i + n` is O(n).
 * Calling *any* non-const member other than `begin()` or `end()` invalidates
 * all iterators! Note that this includes `operator()`, i.e. drawing a random
 * event invalidate all iterators!
 *     
 * There's also a `const_iterator`, which works the same way, but cannot
 * be converted to a `pointer_t`.
 *
 * Event Pointers & EventMovedFunctor:
 * -----------------------------------
 * Event pointers, i.e. `discrete_distribution::pointer`, represents a pointer to
 * an event stored in the distribution. Event pointers aren't plain pointers to
 * instances of `value_type` -- they contain additional information required for
 * updating and removing the event they point to. Any non-const operation other
 * than begin() and end() may cause pointers to existing events to become
 * invalid -- even pointers to events that weren't modified by the operation. If
 * that happens, discrete_distribution invokes the EventMovedFunctor, passing
 * the invalidated event pointer and the new event pointer as arguments. For
 * instances ev_moved_f of EventMovedFunctor, the expression and two values
 * event_invalidated, event_new of type pointer, the expression
 *
 *     ev_moved_f(event_pointer_invalidated, event_pointer_new)
 *
 * must be valid, and should return void (the return value is ignored, though).
 *
 * If no EventMovedFunctor type is specified during discrete_distribution type
 * instantiation, a default no-op functor is used. Such distributions, however,
 * provide *no* way of re-finding an event after a pointer-invalidating
 * operation, other than iterating through the whole distribution. Even then,
 * custom additional data (i.e., a non-void DataType) must be used to identify
 * events, since iteration order is arbitrary.
 *
 * For set() operations, the `EventMovedFunctor` is called if the pointer to the
 * modified event changed as a result of the modification. The return value of
 * set() is the new pointer in that case. For `insert()` and `erase()` operations,
 * the `EventMovedFunctor` is *not* called for the created or removed event, but
 * it will of course be called if the `insert()` or `erase()` operations
 * invalidates pointers to other events.
 *
 * operator() may also trigger EventMovedFunctor invocations. Like set(), if
 * the event that is drawn is also moved, the new pointer is returned.
 *
 * Template Parameter Requirements:
 * ================================
 *
 * SizeType:
 * ---------
 * Mmust be an integral type, and its maximum value limits the maximal
 * number of values that the distribution can contain, i.e. its domain size.
 *
 * WeightType:
 * -----------
 * Must support the following operations, and none of them may throw:
 *
 * ### Default construction (w must afterwards represent the value zero):
 *
 *     WeightType w = WeightType();
 *
 * ### Copy construction:
 *
 *     WeightType w2(w);
 *
 * ### Assignment:
 *
 *     w = w2;
 *
 * ### Zero assignment:
 *
 *     w = 0;
 *
 * ### In-place addition/substraction:
 *
 *     w += w2;
 *     w -= w2;
 *
 * ### Addition/substraction (must return `WeightType`):
 *
 *     w + w2;
 *     w - w2;
 *
 * ### Comparison (must return `bool`):
 *
 *     w <  w_2;
 *     w <= w_2;
 *     w == w_2;
 *     w != w_2;
 *     w >= w_2;
 *     w >  w_2;
 *
 * ### Base-2 logarithm rounded up (must return `std::ptrdiff_t`):
 *
 *     log2ceil(w);
 *
 * ### Power-of-2 for integral exponents of type std::ptrdiff_t:
 *
 *     pow2(exponent, weight_template)
 *
 *
 * DataType:
 * ---------
 * Represents the additional data stored for values. Must be copy constructible,
 * and the copy constructor musn't throw.
 *
 * EventMovedFunctor:
 * ------------------
 * See section on event pointer and `EventMovedFunctor` above.
 *
 * Engine:
 * -------
 * This is not a type parameter of `discrete_distribution`, but rather of its
 * templated `operator()`. `Engine` must support the following operations:
 *
 * ### Uniform weight distribution type instantiation:
 *
 *     typedef dyndist::uniform_distribution<WeightType, Engine> rng_weight_type;
 *     typedef rng_weight_type::result_type rng_weight_result_type;
 *
 * ### Uniform random weight sampling from [0, w):
 *
 *     rng_weight_result_type random_weight = rng_weight_type::sample(w, engine);
 *
 * ### Weight result comparision:
 *
 *     random_weight < w;
 *
 * ### Uniform size distribution type instantiation:
 *
 *     typedef dyndist::uniform_distribution<SizeType, Engine> rng_size_type;
 *     typedef rng_size_type::result_type rng_size_result_type;
 *
 * ### Uniform random size sampling from [0, n):
 *
 *     rng_size_result_type random_size = rng_size_type::sample(n, engine);
 *
 * ### SizeType result conversion:
 *
 *     const SizeType random_size = (SizeType)random;
 *
 * ### Copy Construction & Copy Assignment
 *
 *     rng_weight_result_type r1(random_weight);
 *     rng_size_result_type r2(random_size);
 *     r1 = rng_weight_type::sample(w, engine);
 *     r2 = rng_size_type::sample(n, engine);
 *
 * To reduce the number of calls to the RNG, `Engine` can additionally support
 * the following operations:
 *
 * ### Uniform combined distribution type instantiation:
 *
 *     typedef pow2_expression<WeightType> pow2_expr_t;
 *     typedef std::pair<SizeType, pow2_expr_t> bounds_t;
 *     typedef dyndist::uniform_distribution<bounds_t, Engine> rng_type;
 *     typedef rng_type::result_type rng_pair_result_type;
 *
 * ### Uniform random pair generation within [ 0, n ) x [ 0, 2 ** exp ):
 *
 *     const std::ptrdiff_t exp = ...;
 *     const WeightType two_pow_exp = pow2(exp, weight_template);
 *     const bounds_t bounds = bounds_t(n, pow2_expr_t(exp, two_pow_exp));
 *     rng_pair_result_type random = rng_type::sample(bounds, engine);
 *
 * ### Result comparison and conversion:
 *
 *     const SizeType random_size = (SizeType)random.first;
 *     random.second < w;
 *
 * These operations may throw exceptions.
 *
 * WeightType Considerations:
 * --------------------------
 * Binary operations *don't* necessarily have to be allow arbitrary instances as
 * parameters. Instead, weight type instances may carry a "flavour", and are
 * allowed to produce undefined results if two instances with a different
 * "flavour" are passed to a binary operation. However, no operation may change
 * the "flavour" of its arguments, and the flavour of the result *must* be the
 * same as that of its arguments. This includes the zero-assignment operation,
 * i.e. `w = 0` must assign the value `0` to `w`, but musn't change w's flavour.
 * This implies that *every* flavour of a weight_type must be capable of
 * representing the value zero. Every flavour must also be capable of
 * representing arbitrary integral powers of two, since `pow2()` must also, per
 * the flavour-conservation rule, return an instance with the same flavour as
 * the weight_template argument.
 *
 * The only exception to this flavour-conservation rule is assignment of one
 * weight type to another, i.e. `w = w2`. This operation must change the flavour
 * of `w` to that of `w2`, and assign `w2`'s value. It musn't fail if `w`'s and
 * `w2`'s flavour differ.
 *
 * Random Number Generation Considerations:
 * ----------------------------------------
 * As outlined above, arbitrary random number generators can be used, provided
 * that dyndist::uniform_distribution is specialized for the RNG type (called
 * Engine) and both WeightType and SizeType. Note that the result type of
 * `uniform_distribution<T,Engine>::sample` doesn't have to be `T`! For
 * `T=WeighType`, it sufficies that `r < w` is a valid expressions, and for
 * `T=SizeType` it sufficies that `(SizeType)r` is valid. `T` must additionally
 * be copy constructible and copy assignable.
 *
 * Optional RNG Optimization
 * --------------------------------------------------------
 * In some cases, both a size and a weight are be drawn, and in those cases the
 * upper bound on the weight is always an integral power of two. For integral
 * weight typed, this can be achived with a single draw from the underlying RNG,
 * by drawing a single integer and using its lower bits as the weight, and the
 * upper bits as the size.
 *
 * To facilitate this optimization, discrete_distribution uses
 *
 *     dyndist::uniform_distribution<std::pair<SizeType, pow2_expression<WeightType>,
 *                                           Engine>
 *
 * whenever this optimization is applicable. The default implementation simply
 * obtains the two samples separately by calling
 *
 *     dyndist::uniform_distribution<SizeType>::sample
 *
 * and
 *
 *     dyndist::uniform_distibution<WeightType>::sample,
 *
 * but by providing an explicit specialization of uniform_distribution for this
 * case, one call to the RNG can be avoided as described above.
 */

#ifndef dyndist_discrete_distribution_h
#define dyndist_discrete_distribution_h

#include "dyndist/global.h"

#if DYNDIST_IFSTATS
#include <ostream>
#include <iomanip>
#endif

#include <cstddef>
#include <limits>
#include <algorithm>
#include <cmath>
#include <deque>
#include <vector>

#include "dyndist/utility.h"
#include "dyndist/rangemap.h"

#define DYNDIST_DISCRETE_DISTRIBUTION_FRIEND \
    template<typename,typename,typename,typename> \
    friend class DYNDIST_NAMESPACE::discrete_distribution

DYNDIST_NAMESPACE_BEGIN

/*
 * Forward Declarations
 */
namespace discrete_distribution_details
{
    template<typename WeightType, typename DataType>
    struct value_holder;

    template<typename SizeType, typename WeightType, typename DataType>
    struct level_data;

    template<typename Types>
    struct iterator;
    
    template<typename SizeType, typename WeightType, typename DataType>
    struct f_event_moved_nop;
}

template<
    typename SizeType,
    typename WeightType,
    typename DataType = void,
    typename EventMovedFunctor =
    discrete_distribution_details::f_event_moved_nop<SizeType,WeightType,DataType>
>
class discrete_distribution;


/**
 * \struct discrete_distribution_value
 *
 * Represents events.
 *
 * Events contain a weight of type weight_type, and an optional data member.
 * The value struct is specialized for weight type void. Note that weight is
 * a const member, and values therefore cannot be assigned, but they *can*
 * be copy-constructed.
 */
template<typename WeightType, typename DataType>
struct discrete_distribution_value
{
public:
    typedef WeightType weight_type;
    typedef DataType data_type;
    
    DYNDIST_INLINE
    discrete_distribution_value()
        :weight(), data()
    {}
    
    DYNDIST_INLINE
    discrete_distribution_value
    (const weight_type& _weight, const data_type& _data)
        :weight(_weight), data(_data)
    {}
    
    bool operator==(const discrete_distribution_value& o) const
    { return (weight == o.weight) && (data == o.data); }
    
    const weight_type weight;
    
    data_type data;
    
private:
    DYNDIST_DISCRETE_DISTRIBUTION_FRIEND;
    
    DYNDIST_INLINE
    discrete_distribution_value
    (const discrete_distribution_value& _value, const weight_type& _weight)
        :weight(_weight), data(_value.data)
    {}
};

template<typename WeightType>
struct discrete_distribution_value<WeightType, void>
{
public:
    typedef WeightType weight_type;
    typedef void data_type;
    
    DYNDIST_INLINE
    discrete_distribution_value()
        :weight()
    {}
    
    DYNDIST_INLINE
    discrete_distribution_value(const weight_type& _weight)
        :weight(_weight)
    {}
    
    bool operator==(const discrete_distribution_value& o) const
    { return (weight == o.weight); }
    
    const weight_type weight;
    
private:
    DYNDIST_DISCRETE_DISTRIBUTION_FRIEND;
    
    DYNDIST_INLINE
    discrete_distribution_value
    (const discrete_distribution_value& /* _value */, const weight_type& _weight)
        :weight(_weight)
    {}
};


/**
 * \struct discrete_distribution_pointer
 *
 * Represents a pointer to a value instance.
 *
 * Instead of a raw pointer to the value, pointer instances store a pointer
 * to the level the event is on, as well as the index of the event within
 * the level. Pointers thus remain valid if the std::vector containing an
 * event is resized, but not if events are moved from one level to another.
 * For the latter, discrete_distribution provides a callback interface which
 * notifies clients whenever the pointer to an event changes.
 */
template<typename SizeType, typename WeightType, typename DataType>
struct discrete_distribution_pointer
{
public:
    typedef SizeType size_t;
    typedef WeightType weight_type;
    typedef DataType data_type;
    typedef discrete_distribution_value<weight_type, data_type> value_type;
    
    struct hash
    {
        std::size_t operator()(const discrete_distribution_pointer& p) const;
    };
    
    DYNDIST_INLINE
    discrete_distribution_pointer()
        :level(NULL), index(-size_t(1))
    {}
    
    DYNDIST_INLINE
    operator bool () const
    { return level != NULL; }
    
    DYNDIST_INLINE_FLATTEN
    operator value_type* () const
    { return &(value_type&)as_value_holder(); }
    
    DYNDIST_INLINE_FLATTEN
    value_type* operator->() const
    { return &(value_type&)as_value_holder(); }
    
    DYNDIST_INLINE_FLATTEN
    value_type& operator*() const
    { return as_value_holder(); }
    
    DYNDIST_INLINE
    bool operator< (const discrete_distribution_pointer& o) const
    { return (level != o.level) ? (level < o.level) : (index <  o.index); }
    
    DYNDIST_INLINE
    bool operator<=(const discrete_distribution_pointer& o) const
    { return (level != o.level) ? (level < o.level) : (index <= o.index); }
    
    DYNDIST_INLINE
    bool operator>=(const discrete_distribution_pointer& o) const
    { return (level != o.level) ? (level > o.level) : (index >= o.index); }
    
    DYNDIST_INLINE
    bool operator> (const discrete_distribution_pointer& o) const
    { return (level != o.level) ? (level > o.level) : (index >  o.index); }
    
    DYNDIST_INLINE
    bool operator==(const discrete_distribution_pointer& o) const
    { return (level == o.level) && (index == o.index); }
    
    DYNDIST_INLINE
    bool operator!=(const discrete_distribution_pointer& o) const
    { return (level != o.level) || (index != o.index); }
    
private:
    DYNDIST_DISCRETE_DISTRIBUTION_FRIEND;
    
    template<typename>
    friend struct discrete_distribution_details::iterator;
    
    typedef discrete_distribution_details::level_data
    <size_t, weight_type, data_type>
    level_data_type;
    
    typedef std::pair<std::ptrdiff_t,  level_data_type>
    levels_value_type;
    
    typedef discrete_distribution_details::value_holder<weight_type, data_type>
    level_events_value_type;
    
    DYNDIST_INLINE
    discrete_distribution_pointer(levels_value_type& _level, value_type& _event)
        :level(&_level)
        ,index((size_t)(&_event - &(value_type&)_level.second.events.front()))
    {
        DYNDIST_ASSERT(index < _level.second.events.size());
        DYNDIST_ASSERT(*this == &_event);
    }
    
    DYNDIST_INLINE
    level_events_value_type& as_value_holder() const
    { return level->second.events[(std::size_t)index]; }
    
    levels_value_type* level;
    size_t index;
};

template<typename S, typename W, typename D>
DYNDIST_INLINE
std::size_t
discrete_distribution_pointer<S,W,D>::hash::operator()
(const discrete_distribution_pointer& p) const
{
    const std::size_t N = sizeof(std::size_t)*8;
    const std::size_t l = reinterpret_cast<char*>(p.level) - (char*)NULL;
    return (l << (N/2)) ^ (l >> (N/2)) ^ (std::size_t)p.index;
}


/**
 * Contains helper types for discrete_distribution
 */
namespace discrete_distribution_details
{
    /**
     * \typedef level_t
     *
     * Represents weight levels, i.e. ceil(log2(weight)).
     */
    typedef std::ptrdiff_t level_t;
    

    /**
     * \struct value_holder
     *
     * Used to store values in a std::vector.
     *
     * The standard allows std::vector to demand a functional assignment
     * operator, even if the only vector operations used are push_back() and
     * pop_back(). This is e.g. the case for libg++, at least in C++03 mode.
     * Since weight is a const member of value, this prevents value from
     * directly being used as a std::vector element type. The value_holder works
     * around that by providing an assignment operator which destructs the old
     * value, and constructs the new value in place. Doing so does not defeat
     * the purpose of making weight a const member, because
     * discrete_distribution's API is carefull to not leak any references to
     * value_holder instances, but only to the value instances they contain.
     */
    template<typename WeightType, typename DataType>
    struct value_holder
    {
        typedef WeightType weight_type;
        
        typedef DataType data_type;
        
        typedef discrete_distribution_value<weight_type, data_type> value_type;
        
        DYNDIST_INLINE
        operator value_type& ()
        { return m_value; }

        DYNDIST_INLINE
        operator const value_type& () const
        { return m_value; }

        DYNDIST_INLINE
        value_type& operator=(const value_holder& value)
        { m_value.~value_type(); new(&m_value) value_type(value); return *this; }
        
    private:
        DYNDIST_DISCRETE_DISTRIBUTION_FRIEND;

        DYNDIST_INLINE
        value_holder(const value_type& value)
            :m_value(value)
        {}
        
        value_type m_value;
    };
    
    
    template<typename DistributionType>
    struct iterator;

    
    /**
     * \struct level_data
     *
     * Represents a weight level, and physically stores the events beloning
     * to such a level.
     */
    template<typename SizeType, typename WeightType, typename DataType>
    struct level_data
    {
        typedef SizeType size_t;
        
        typedef WeightType weight_type;
        
        typedef DataType data_type;
        
        typedef std::pair<std::ptrdiff_t, level_data> levels_value_type;
        
        typedef rangemap<level_data> levels_type;
        
        typedef discrete_distribution_value<weight_type, data_type> value_type;

        typedef value_holder<weight_type, data_type> value_holder_type;

        typedef std::vector<value_holder_type> events_type;

        typedef discrete_distribution_value<weight_type, void*>
        meta_value_type;
        
        typedef discrete_distribution_pointer<size_t, weight_type, void*>
        meta_pointer;

        /* Event-moved callback for the meta distribution.
         * Updates the level's link_meta field if the meta event moves
         */
        struct f_meta_moved
        {
            void operator()(meta_pointer event, meta_pointer event_new);
        };

        /* Functed used to initialize a level */
        template<typename DistributionType>
        struct f_level_init
        {
            f_level_init(DistributionType& dist)
                :m_dist(&dist)
            {}
            
            void operator()(levels_value_type& level, const bool front);
            
        private:
            DistributionType* m_dist;
        };

        DYNDIST_INLINE
        level_data()
            :link_nonzero(-std::size_t(1)), link_updated(-std::size_t(1))
            ,link_meta()
            ,weight_lower(), weight_upper(), total_weight()
            ,events()
            ,tiny(false), pure(true)
        {}
        
        /** Position of the level within distribution's nonzero-level list, or -1 */
        std::size_t link_nonzero;

        /** Position of the level within distribution's updated-level list, or -1 */
        std::size_t link_updated;
        
        /** Pointer to the level's event in the meta distribution */
        meta_pointer link_meta;

        /**
         * Lower and upper bound on the level's event.
         * All events on a level always satisfy weight <= weight_upper.
         * All events on a *pure* level also satisfy weight > weight_lower.
         */
        weight_type weight_lower;
        weight_type weight_upper;
        
        /** Sum of weights of the level's events */
        weight_type total_weight;
        
        /** Container holding the level's events */
        events_type events;
        
        /** True if the level is currently the level with the smallest lower bound */
        bool tiny;
         
        /** True if weight > weight_lower holds for all events on the level */
        bool pure;
    };
    
    template<typename S, typename W, typename D>
    DYNDIST_INLINE
    void
    level_data<S,W,D>::f_meta_moved::operator()
#if DYNDIST_IFASSERT
    (meta_pointer event, meta_pointer event_new)
#else
    (meta_pointer, meta_pointer event_new)
#endif
    {
        DYNDIST_ASSERT(event->data == event_new->data);
        
        const meta_value_type& event_value = *event_new;
        levels_value_type& level =
        *reinterpret_cast<levels_value_type*>(event_value.data);
        level_data& ld = level.second;
        DYNDIST_ASSERT(ld.link_meta == event);
        ld.link_meta = event_new;
    }
    
    template<typename S, typename W, typename D>
    template<typename D2>
    DYNDIST_INLINE
    void
    level_data<S,W,D>::f_level_init<D2>::operator()
    (levels_value_type& level, const bool)
    {
        const level_t lk = level.first;
        level_data& ld = level.second;
        
        /* Set level's weight and weight bounds */
        ld.weight_lower = pow2(lk - 1, m_dist->m_zero_weight);
        ld.weight_upper = pow2(lk, m_dist->m_zero_weight);
        ld.total_weight = m_dist->m_zero_weight;
        
        /* Add corresponding meta event and link with level */
        DYNDIST_ASSERT(level.second.link_meta == meta_pointer());
        
        if (m_dist->m_meta)
            level.second.link_meta =
            m_dist->m_meta->insert(meta_value_type(m_dist->m_zero_weight,
                                                   &level));
    }
    
    
    /**
     * \struct iterator
     *
     * A random access iterator which allows iteration through all the events
     * stored within a discrete_distribution.
     *
     * Iterators are somewhat similar to pointers in that they reference a
     * particular event within a discrete_distribution, but have much larger
     * overhead than these. Pointer, however, provide no means of moving from
     * one event to the next or previous one.
     *
     * Calling *any* non-const member of discret_distribution potentially
     * invalidates *all* iterators. Note that this includes drawing elements
     * from a discrete_distribution -- operator() is non-cost!
     *
     * The iterator produced by end() has the last level as its current level,
     * and m_event_i is the end iterator of that level's events vector. Note
     * that this doesn't imply that m_level_i of such an iterator is
     * dereferencable -- the distribution may *only* contain a zero level, so
     * m_level_i might equal m_levels_end.
     *
     * The iterator produced by begin() has the first non-empty level as its
     * current level, if such a level exists, and equals end() otherwise. In
     * the first case, m_event_i points to the first element of the first
     * non-empty level.
     */
    template<typename Types>
    struct iterator
    {
        typedef typename Types::distribution_type distribution_type;

        typedef typename Types::value_type value_type;

        typedef typename Types::pointer pointer;
        
        typedef typename distribution_type::size_t size_t;

    private:
        DYNDIST_DISCRETE_DISTRIBUTION_FRIEND;

        typedef typename Types::levels_iterator levels_iterator;
        typedef typename Types::events_iterator events_iterator;
        typedef typename Types::levels_value_type levels_value_type;
        
        static
        iterator begin(distribution_type& distribution);

        static
        iterator end(distribution_type& distribution);

        DYNDIST_INLINE_FLATTEN
        iterator
        (const levels_iterator& level_i, levels_value_type& level,
         const events_iterator& event_i, distribution_type& distribution)
            :m_distribution(&distribution)
            ,m_level_i(level_i)
            ,m_level(&level)
            ,m_event_i(event_i)
        {
            DYNDIST_ASSERT(m_level_i >= m_distribution->m_levels_nonzero.begin());
            DYNDIST_ASSERT(m_level_i <= m_distribution->m_levels_nonzero.end());
            DYNDIST_ASSERT(m_level ==
                         ((m_level_i != m_distribution->m_levels_nonzero.end())
                          ? *m_level_i
                          : &m_distribution->m_level_zero));
            DYNDIST_ASSERT(m_event_i >= m_level->second.events.begin());
            DYNDIST_ASSERT(m_event_i <= m_level->second.events.end());
        }

    public:
        DYNDIST_INLINE_FLATTEN
        iterator()
            :m_distribution(NULL)
            ,m_level_i(), m_level(NULL), m_event_i()
        {}
        
        DYNDIST_INLINE_FLATTEN
        operator pointer () const
        { return pointer(*m_level, *m_event_i); }
        
        DYNDIST_INLINE_FLATTEN
        value_type& operator*() const
        { return *m_event_i; }
        
        DYNDIST_INLINE_FLATTEN
        value_type* operator->() const
        { return &(value_type&)*m_event_i; }
        
        DYNDIST_INLINE
        value_type& operator[](std::ptrdiff_t offset) const
        { return *(*this + offset); }
        
        DYNDIST_INLINE
        iterator& operator++()
        { ++m_event_i; adjust_fwd(); return *this; }

        DYNDIST_INLINE
        iterator& operator+=(std::ptrdiff_t d)
        { (d > 0) ? fwd((size_t)d) : rev((size_t)-d); return *this; }

        DYNDIST_INLINE
        iterator& operator--()
        { adjust_rev(); --m_event_i; return *this; }

        DYNDIST_INLINE
        iterator& operator-=(std::ptrdiff_t d)
        { (d > 0) ? rev((size_t)d) : fwd((size_t)-d); return *this; }

        DYNDIST_INLINE
        iterator operator++(int) const
        { iterator i = *this; ++(*this); return i; }

        DYNDIST_INLINE
        iterator operator+(std::ptrdiff_t d) const
        { iterator i = *this; i += d; return i; }

        DYNDIST_INLINE
        iterator operator--(int) const
        { iterator i = *this; --(*this); return i; }

        DYNDIST_INLINE
        iterator operator-(std::ptrdiff_t d) const
        { iterator i = *this; i -= d; return i; }

        DYNDIST_INLINE
        std::ptrdiff_t operator-(const iterator& o) const
        { return ((*this < o)
                  ? -(std::ptrdiff_t)dist_fwd(o)
                  :  (std::ptrdiff_t)o.dist_fwd(*this)); }
        
        DYNDIST_INLINE
        bool operator< (const iterator& o) const
        { return ((m_level_i < o.m_level_i)
                  ? true
                  : ((m_level_i == o.m_level_i)
                     ? m_event_i < o.m_event_i
                     : false)); }

        DYNDIST_INLINE
        bool operator<= (const iterator& o) const
        { return ((m_level_i < o.m_level_i)
                  ? true
                  : ((m_level_i == o.m_level_i)
                     ? m_event_i <= o.m_event_i
                     : false)); }

        DYNDIST_INLINE_FLATTEN
        bool operator==(const iterator& o) const
        { return (m_level == o.m_level) && (m_event_i == o.m_event_i); }

        DYNDIST_INLINE_FLATTEN
        bool operator!=(const iterator& o) const
        { return !(*this == o); }

        DYNDIST_INLINE_FLATTEN
        bool operator>=(const iterator& o) const
        { return !(*this < o); }

        DYNDIST_INLINE_FLATTEN
        bool operator> (const iterator& o) const
        { return !(*this <= o); }

    private:
        /** Move delta entries forward  */
        void fwd(size_t delta);

        /** Move delta entries backward */
        void rev(size_t delta);

        /** Measure distance towards b, which must be larger */
        size_t dist_fwd(const iterator& b) const;

        /** Move from end of level to beginning of next level */
        void adjust_fwd();
        
        /** Move from beginning of level to end of previous level */
        void adjust_rev();
        
        /** Move to first entry on the next level */
        void level_fwd();
        
        /** Move to end of previous level */
        void level_rev();

        /** The distribution */
        distribution_type* m_distribution;

        /**
         * The current level iterator.
         *
         * m_level == *m_level_i if m_level_i is dereferencable,
         * m_level == &m_distribution->m_level_zero otherwise.
         */
        levels_iterator m_level_i;

        /** The current level */
        levels_value_type* m_level;
        
        /** The current event */
        events_iterator m_event_i;
    };

    template<typename D>
    DYNDIST_FLATTEN
    iterator<D>
    iterator<D>::begin(distribution_type& distribution)
    {
        if (distribution.m_levels_nonzero.empty())
            /* Return iterator positioned at first zero-level event */
            return iterator(distribution.m_levels_nonzero.end(),
                            distribution.m_level_zero,
                            distribution.m_level_zero.second.events.begin(),
                            distribution);
        else {
            /* Return iterator positioned at first nonzero-level event */
            const levels_iterator first_i = distribution.m_levels_nonzero.begin();
            levels_value_type& first = **first_i;
            return iterator(first_i, first, first.second.events.begin(),
                            distribution);
        }
    }

    template<typename D>
    DYNDIST_FLATTEN
    iterator<D>
    iterator<D>::end(distribution_type& distribution)
    {
        /* Return iterator positioned at zero-level end */
        return iterator(distribution.m_levels_nonzero.end(),
                        distribution.m_level_zero,
                        distribution.m_level_zero.second.events.end(),
                        distribution);
    }

    template<typename D>
    DYNDIST_FLATTEN
    void
    iterator<D>::fwd(size_t delta)
    {
        /* Skip levels until delta becomes small than the current level's size.
         * The test for delta > 0 is necessary due to the zero level, which
         * may be empty. If we reach the zero level at the same time that delta
         * reaches zero, we'd continue without the additional test.
         */
        for(size_t level_remaining = (size_t)(m_level->second.events.end() - m_event_i);
            (delta > 0) && (delta >= level_remaining);
            level_remaining = (size_t)(m_level->second.events.end() - m_event_i))
        {
            /* Invariants
             * level_remaining > 0 could only be violated if the zero-level is
             * empty. But then, delta should have reached zero by the time we
             * reach the zero level, since we're otherwise trying to forward
             * past the end iterator, which is illegal.
             */
            DYNDIST_ASSERT(delta > 0);
            DYNDIST_ASSERT(level_remaining > 0);
            
            /* Skip the remaining events on the current level */
            delta -= level_remaining;
            level_fwd();
        }
        
        /* Skip remaining events within the current level */
        m_event_i += (std::ptrdiff_t)delta;
    }

    template<typename D>
    DYNDIST_FLATTEN
    void
    iterator<D>::rev(size_t delta)
    {
        for(size_t level_remaining = (size_t)(m_event_i - m_level->second.events.begin());
            delta > level_remaining;
            level_remaining = (size_t)(m_event_i - m_level->second.events.begin()))
        {
            /* Invariants */
            DYNDIST_ASSERT(delta > 0);
            
            /* Skip the remaining events on the current level */
            delta -= level_remaining;
            level_rev();
        }
        
        /* Skip remaining events within the current level */
        m_event_i -= (std::ptrdiff_t)delta;
    }

    template<typename D>
    DYNDIST_FLATTEN
    typename iterator<D>::size_t
    iterator<D>::dist_fwd(const iterator& b) const
    {
        DYNDIST_ASSERT(*this <= b);
        size_t d=0;
        
        /* Account for elements on the levels up to b's level */
        iterator a = *this;
        while(a.m_level_i != b.m_level_i) {
            d += (size_t)(a.m_level->second.events.end() - a.m_event_i);
            a.level_fwd();
        }
        
        /* Account for elements on b's level */
        d += (size_t)(b.m_event_i - a.m_event_i);
        
        return d;
    }

    template<typename D>
    DYNDIST_INLINE
    void
    iterator<D>::adjust_fwd()
    {
        /* Move from end of a non-zero level to beginning of next level */
        if ((m_event_i == m_level->second.events.end()) &&
            (m_level != &m_distribution->m_level_zero))
            level_fwd();
    }

    template<typename D>
    DYNDIST_INLINE
    void
    iterator<D>::adjust_rev()
    {
        /* Move from beginning of a level to end of the previous non-zero level */
        if (m_event_i == m_level->second.events.begin())
            level_rev();
    }
    
    template<typename D>
    DYNDIST_FLATTEN
    void
    iterator<D>::level_fwd()
    {
        /* Move to the next level */
        DYNDIST_ASSERT(m_level_i != m_distribution->m_levels_nonzero.end());
        ++m_level_i;
        
        /* If we reached the end of the non-zero levels, use zero level */
        if (m_level_i != m_distribution->m_levels_nonzero.end())
            m_level = *m_level_i;
        else
            m_level = &m_distribution->m_level_zero;
        m_event_i = m_level->second.events.begin();
    }

    template<typename D>
    DYNDIST_FLATTEN
    void
    iterator<D>::level_rev()
    {
        /* Move to the previous level */
        DYNDIST_ASSERT(m_level_i != m_distribution->m_levels_nonzero.begin());
        --m_level_i;
        m_level = *m_level_i;
        m_event_i = m_level->second.events.end();
    }
    
    
    /**
     * \struct iterator_types
     *
     * Types for discrete_distribution non-const iterators
     */
    template<typename DistributionType>
    struct iterator_types {
        typedef DistributionType distribution_type;
        typedef typename distribution_type::value_type value_type;
        typedef typename distribution_type::pointer pointer;
        typedef typename distribution_type::levels_collection_type::iterator
        levels_iterator;
        typedef typename distribution_type::level_events_type::iterator
        events_iterator;
        typedef typename distribution_type::levels_value_type
        levels_value_type;
    };
    
    
    /**
     * \struct const_iterator_types
     *
     * Types for discrete_distribution non-const iterators
     */
    template<typename DistributionType>
    struct const_iterator_types {
        typedef const DistributionType distribution_type;
        typedef const typename distribution_type::value_type value_type;
        typedef typename distribution_type::pointer pointer;
        typedef typename distribution_type::levels_collection_type::const_iterator
        levels_iterator;
        typedef typename distribution_type::level_events_type::const_iterator
        events_iterator;
        typedef const typename distribution_type::levels_value_type
        levels_value_type;
    };
    

    /**
     * \struct f_event_moved_nop
     *
     * Default functor for the event-moved callback parameters. Does nothing.
     */
    template<typename SizeType, typename WeightType, typename DataType>
    struct f_event_moved_nop
    {
        typedef discrete_distribution_pointer<SizeType, WeightType, DataType>
        pointer;
        
        DYNDIST_INLINE
        void operator()(const pointer&, const pointer&)
        {}
    };
}


/*
 * discrete_distribution
 *
 * See documentation at the top
 */
template<typename SizeType, typename WeightType, typename DataType,
         typename EventMovedFunctor>
class discrete_distribution
{
    /*
     * Template Parameters & Public Types
     */
public:
    typedef SizeType size_t;
    typedef WeightType weight_type;
    typedef DataType data_type;
    typedef EventMovedFunctor event_moved_functor;

    /*
     * Value & Pointer Types
     */
public:
    typedef discrete_distribution_value<weight_type, data_type>
    value_type;

    typedef discrete_distribution_pointer<size_t, weight_type, data_type>
    pointer;

    friend struct discrete_distribution_details::iterator_types
    <discrete_distribution>;
    friend struct discrete_distribution_details::iterator
    <discrete_distribution_details::iterator_types<discrete_distribution> >;
    
    typedef discrete_distribution_details::iterator
    <discrete_distribution_details::iterator_types<discrete_distribution> >
    iterator;

    friend struct discrete_distribution_details::const_iterator_types
    <discrete_distribution>;
    friend struct discrete_distribution_details::iterator
    <discrete_distribution_details::const_iterator_types<discrete_distribution> >;

    typedef discrete_distribution_details::iterator
    <discrete_distribution_details::const_iterator_types<discrete_distribution> >
    const_iterator;
    
    /*
     * Meta-Distribution & Levels
     */
private:
    typedef std::ptrdiff_t level_type;
    
    /* Level Data */
    typedef discrete_distribution_details::level_data
    <size_t, weight_type, data_type>
    level_data_type;
    
    typedef typename level_data_type::levels_value_type levels_value_type;
    typedef typename level_data_type::levels_type levels_type;
    
    /* Event Data */
    typedef typename level_data_type::events_type level_events_type;
    typedef typename level_events_type::value_type level_events_value_type;
    
    /* Level Collections (used for nonzero- and updated-levels collection) */
    typedef std::vector<levels_value_type*> levels_collection_type;
    
    /* Meta Distribution
     *
     * Uses void* instead of levels_value_type* as the DataType to avoid an
     * infinite template depth. If we used the correct type, each meta-level
     * would require a distinct level_data_type, because starting from the
     * n-th level levels_value_type, the expression
     *   *levels_value_type::events_type::element_type::data_type
     * yields the (n-1)-th level levels_value_type, or the distribution's
     * actual data_type if n=0. Some kind of type erasure is thus necessary.
     */
    typedef discrete_distribution
    <size_t, weight_type, void*, typename level_data_type::f_meta_moved>
    meta_type;
    
    typedef typename meta_type::value_type meta_value_type;
    
    typedef typename meta_type::pointer meta_pointer;
    
    /* Level initialization functor */
    template<typename> friend struct level_data_type::f_level_init;
    
    typedef typename level_data_type::template f_level_init
    <discrete_distribution>
    f_level_init;
    
    /*
     * Construction & Destruction
     */
public:
    discrete_distribution(const event_moved_functor& event_moved =
                          event_moved_functor());

    ~discrete_distribution();

private:
    /** discrete_distribution instances are non-copyable */
    discrete_distribution(const discrete_distribution& other);

    /** discrete_distribution instances are non-assignable */
    void operator= (const discrete_distribution& other);
    
    /*
     * Observers
     */
public:
    static const size_t max_size;
    
    /** Returns true if there are no events with non-zero weight */
    DYNDIST_INLINE
    bool empty() const
    { return m_nonzero_count == 0; }
    
    /** Returns the number of events */
    DYNDIST_INLINE
    size_t size() const
    { return m_total_count; }

    /** Returns the number of events with non-zero weight */
    DYNDIST_INLINE
    size_t size_nonzero() const
    { return m_nonzero_count; }

    /** Returns the sum of the event's weights */
    DYNDIST_INLINE
    const weight_type& weight() const
    { return m_total_weight; }

    /** Returns an iterator to the first event, or end() if size() == 0 */
    DYNDIST_INLINE
    iterator begin()
    { return iterator::begin(*this); }

    /** Returns an iterator to the first event, or end() if size() == 0 */
    DYNDIST_INLINE
    const_iterator begin() const
    { return const_iterator::begin(*this); }

    /** Returns the end iterator */
    DYNDIST_INLINE
    iterator end()
    { return iterator::end(*this); }

    /** Returns the end iterator */
    DYNDIST_INLINE
    const_iterator end() const
    { return const_iterator::end(*this); }

    /*
     * Mutators
     */
public:
    /**
     * Erases all values
     *
     * Invalidates iterators and pointers.
     */
    void clear();

    /**
     * Inserts a new value
     *
     * \param value the value to insert
     * \return a pointer to the inserted value
     *
     * Invalidates iterators and pointers.
     */
    pointer insert(const value_type& value);
    
    /**
     * Erases an event
     *
     * Invalidates iterators and pointers.
     */
    void erase(pointer event);
    
    /**
     * Erases an event
     *
     * \param pointer to the event to update
     * \param value new eight
     * \return pointer to the updated event
     *
     * Invalidates iterators and pointers.
     */
    pointer set(pointer event, const weight_type& value);

    /**
     * \struct batch_update
     *
     * RAII-style batched_update proxy.
     */
    struct batch_update
    {
        /**
         * Start a batch of updates
         *
         * \param distribution the `discrete_distribution` that is updated
         * \param total_count_final final number of values in the distribution
         * \param nonzero_count_final final number of values with non-zero weight
         * \param weight_final final total weight of the distribution
         *
         * For a particular `discrete_distribution` instance, only one
         * `batch_update` instance may exist at any point in time.
         */
        batch_update(discrete_distribution& distribution,
                     const size_t total_count_final,
                     const size_t nonzero_count_final,
                     const weight_type& weight_final);
        
        /**
         * Completes the batch of updates.
         */
        ~batch_update();
        
        /**
         * Inserts a new value
         *
         * \param value the value to insert
         * \return a pointer to the inserted value
         *
         * Invalidates iterators and pointers.
         */
        pointer insert(const value_type& value);
        
        /**
         * Erases an event
         *
         * Invalidates iterators and pointers.
         */
        void erase(pointer event);
        
        /**
         * Erases an event
         *
         * \param pointer to the event to update
         * \param value new eight
         * \return pointer to the updated event
         *
         * Invalidates iterators and pointers.
         */
        pointer set(pointer event, const weight_type& value);
        
    private:
        /* Prevent Copying & Assignment */
        batch_update(const batch_update&);
        batch_update& operator=(const batch_update&);
                     
        discrete_distribution& m_distribution;
        bool m_utc;
#if DYNDIST_IFASSERT
        size_t m_total_count;
        size_t m_nonzero_count;
        weight_type m_total_weight;
#endif
    };
    
    /** Returns the current force-meta setting */
    DYNDIST_FLATTEN
    char force_meta() const
    { return m_force_meta; }
    
    /** Pass 1 to always use the meta distribution, -1 to neve use it, 0 for auto */
    DYNDIST_FLATTEN
    void force_meta(char force_meta)
    { m_force_meta = force_meta; }

    /*
     * Sampling
     */
public:
    /**
     * Select a random value, each value with a probability proportional to
     * its weight.
     *
     * \param engine the RNG to use
     * \return a pointer the selected value.
     *
     * Invalidates iterators and pointers.
     */
    template<typename Engine>
    pointer operator()
    (Engine& engine);

    /**
     * Select a random value, each value with a probability proportional to
     * its weight.
     *
     * \param engine the RNG to use
     * \param a reference to a uniformly distributed value in `[ 0, weight() )`.
     *        Upon return, `x` contains a uniformy distributed value in
     *        `[ 0, s->weight )` where `s` is the returned pointer.
     * \return a pointer the selected value.
     *
     * Invalidates iterators and pointers.
     */
    template<typename Engine>
    pointer operator()
    (Engine& engine,
     typename uniform_distribution<weight_type, Engine>::result_type& x);

    
    /*
     * Statistics
     */
public:
#if DYNDIST_IFSTATS
    void print_statistics(std::ostream& dst, std::size_t meta_depth = 0) const;
    
    std::size_t stats_draws() const
    { return m_draws; }
    
    std::size_t stats_draws_level_meta() const
    { return m_draws_level_meta; }
    
    std::size_t stats_draws_level_scan() const
    { return m_draws_level_scan; }
    
    std::size_t stats_draws_level_scan_iterations() const
    { return m_draws_level_scan_iterations; }
    
    std::size_t stats_draws_single() const
    { return m_draws_single; }
    
    std::size_t stats_draws_acceptor() const
    { return m_draws_acceptor; }
    
    std::size_t stats_draws_acceptor_iterations() const
    { return m_draws_acceptor_iterations; }

    std::size_t stats_draws_scan() const
    { return m_draws_scan; }

    std::size_t stats_draws_scan_iterations() const
    { return m_draws_scan_iterations; }

    std::size_t stats_draws_purify() const
    { return m_draws_purify; }
#endif
    
    /*
     * Verification
     */
#if DYNDIST_IFVERIFY
    void verify_consistency() const;
    
    void verify_meta_consistency();
#endif

    /*
     * Private Members
     */
private:
    /**
     * Draws a level by scanning m_levels_nonzero.
     * X is some uniform_distribution<weight_type, Engine>::result_type.
     */
    template<typename X>
    levels_value_type& draw_level_scan(X& x);

    /**
     * Draws a level using the meta distribution.
     */
    template<typename Engine>
    levels_value_type&
    draw_level_meta
    (Engine& engine,
     typename uniform_distribution<weight_type, Engine>::result_type& x);

    /**
     * Draws an event from a level with the acceptor method.
     */
    template<typename Engine>
    pointer
    draw_acceptor
    (levels_value_type& level_data, Engine& engine,
     typename uniform_distribution<weight_type, Engine>::result_type& x);

    /**
     * Draws an event from a level by scanning the level's events.
     * X is some uniform_distribution<weight_type, Engine>::result_type.
     */
    template<typename X>
    pointer draw_scan(levels_value_type& level_data, X& x);

    /**
     * Draws an event from a level by scanning the level's events. Additionally
     * moves all the events which don't exceed the level's lower bound off the
     * level, which allows draw_acceptor to be used the *next time* an event
     * is to be selected from this level.
     * X is some uniform_distribution<weight_type, Engine>::result_type.
     */
    template<typename X>
    pointer draw_purify(levels_value_type& level_data, X& x);

    /**
     * Inserts an event as part of a batched update.
     */
    pointer batch_insert(const value_type& value);
    
    /**
     * Updates an event's weight as part of a batched update.
     */
    pointer batch_set(pointer event, const weight_type& value);

    /**
     * Removes an event as part of a batched update.
     */
    void batch_erase(pointer event);

    /**
     * Updates an event's weight.
     * update_totals() must already have run.
     */
    pointer update_weight(pointer event, value_type& event_value,
                          const weight_type& weight,
                          bool weight_zero, bool weight_was_zero);

    /**
     * Updates an event's weight from zero to a non-zero value.
     * update_totals() must already have run.
     */
    pointer update_weight_zero_nonzero(pointer event,
                                       const weight_type& weight);

    /**
     * Updates an event's weight from a non-zero to a non-zero value.
     * update_totals() must already have run.
     */
    pointer update_weight_nonzero_nonzero(pointer event,
                                          const weight_type& weight);

    /**
     * Updates an event's weight from a non-zero value to zero.
     * update_totals() must already have run.
     */
    pointer update_weight_nonzero_zero(pointer event_ptr);
    
    /**
     * Adds an event to the specified level, which musn't be the zero level.
     */
    pointer place(const value_type& event, levels_value_type& level);

    /**
     * Add an event to the zero level
     */
    pointer place_zero(const value_type& event);

    /**
     * Removes an event from its current level, which musn't be the zero level.
     */
    void unplace(pointer event);

    /*
     * Removes an event from the zero level.
     */
    void unplace_zero(pointer event);
    
    /**
     * Determine the level on which the specified weight belongs.
     */
    level_type weight_level(const weight_type& weight);
    
    /**
     * Setup things upon seeing the first weight_type instance.
     */
    void initialize_weight_type(const weight_type& weight);
    
    /**
     * Update the distributions total weight and number of non-zero events.
     */
    bool update_totals(size_t nonzero_count, const weight_type& weight);

    /**
     * update_totals()'s slow path.
     * This is a separate function to allow the fast path to be inlined.
     */
    void update_totals_begin(bool weight_within_slack,
                             bool count_within_slack);

    /**
     * Cleanup levels after update_totals().
     * Must be called whenever update_totals() returns true, and before
     * events are drawn from the distribution. update_totals() may be called
     * multiple times, and multiple events may be inserted, updated or removed,
     * before update_totals_complete() is finally called. So long, that is,
     * as that call happens before any call to operator(), i.e. before events
     * are drawn from the distribution.
     */
    void update_totals_complete();

    /**
     * Removes the specified level from the non-zero levels collection.
     */
    void levels_nonzero_remove(levels_value_type& level);

    /**
     * Adds the specified level to the updated-levels collection.
     */
    void levels_updated_add(levels_value_type& level);

    /**
     * Removes the specified level from the updated-levels collection.
     */
    void levels_updated_remove(levels_value_type& level);

    /**
     * Updates (and creates, if necessary) the meta distribution.
     */
    void update_meta();
    
    /**
     * Private Data
     */
private:
    /**
     * Threshold for switching from draw_level_scan to draw_level_meta
     */
    static const std::size_t meta_threshold = 16;

    /**
     * Threshold for switching from draw_scan to draw_acceptor
     */
    static const std::size_t scan_threshold = 8;

    /**
     * The level index corresponding to the zero level.
     * Note that the zero level isn't part of m_levels. This constants only
     * serves as a convenient initial value for m_level_tiny and m_level_total.
     */
    static const level_type zero_level;
    
    /**
     * The functor that is called whenever an event moves, i.e. if the pointer
     * to an event is changes.
     */
    event_moved_functor m_event_moved;
    
    /**
     * The weight_type representation of "0".
     */
    weight_type m_zero_weight;
    
    /**
     * The total number of events
     */
    size_t m_total_count;

    /**
     * The number of events with non-zero weight
     */
    size_t m_nonzero_count;
    
    /**
     * The total weight, i.e. the sum of all the event's weights.
     */
    weight_type m_total_weight;
    
    /**
     * The bounds within m_nonzero_count may change without requiring levels
     * to be created or collapsed.
     */
    size_t m_nonzero_count_lower;
    size_t m_nonzero_count_upper;
    
    /**
     * The bounds within m_total_weight may change without requiring levels
     * to be created or collapsed.
     */
    weight_type m_total_weight_lower;
    weight_type m_total_weight_upper;
    
    /**
     * Whether or not some call to update_totals() returned true since the
     * last call to update_totals_complete()
     */
#if DYNDIST_IFASSERT
    bool m_update_totals_needs_completion;
#endif
    
    /**
     * Whether there is a `batch_update` operation in progress
     */
#if DYNDIST_IFASSERT
    bool m_batch_update_in_progress;
#endif
    
    /**
     * The rangemap mapping level indices to level_data_type instances.
     */
    levels_type m_levels;
    
    /**
     * The level data for the zero level.
     */
    levels_value_type m_level_zero;
    
    /**
     * The collection of levels with non-zero weight
     */
    levels_collection_type m_levels_nonzero;
    
    /**
     * The collection of levels whose weight doesn't match their meta events
     */
    levels_collection_type m_levels_updated;
    
    level_type m_level_total;
    
    level_type m_level_tiny;
    
    char m_force_meta;
    
    /**
     * The meta distribution
     */
    meta_type* m_meta;
    
    /**
     * Statistics
     */
private:
#if DYNDIST_IFSTATS
    mutable std::size_t m_draws;
    mutable std::size_t m_draws_level_meta;
    mutable std::size_t m_draws_level_scan;
    mutable std::size_t m_draws_level_scan_iterations;
    mutable std::size_t m_draws_single;
    mutable std::size_t m_draws_acceptor;
    mutable std::size_t m_draws_acceptor_iterations;
    mutable std::size_t m_draws_scan;
    mutable std::size_t m_draws_scan_iterations;
    mutable std::size_t m_draws_purify;
#endif
};

template<typename S,typename W,typename D,typename C>
const typename discrete_distribution<S,W,D,C>::size_t
discrete_distribution<S,W,D,C>::max_size =
std::numeric_limits<typename discrete_distribution<S,W,D,C>::size_t>::max() - 1;

template<typename S,typename W,typename D,typename C>
const typename discrete_distribution<S,W,D,C>::level_type
discrete_distribution<S,W,D,C>::zero_level =
std::numeric_limits<typename discrete_distribution<S,W,D,C>::level_type>::min();

template<typename S,typename W,typename D,typename C>
const std::size_t
discrete_distribution<S,W,D,C>::meta_threshold;

template<typename S,typename W,typename D,typename C>
const std::size_t
discrete_distribution<S,W,D,C>::scan_threshold;

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
discrete_distribution<S,W,D,C>::discrete_distribution
(const event_moved_functor& event_moved)
    :m_event_moved(event_moved)
    ,m_zero_weight()
    ,m_total_count(0), m_nonzero_count(0), m_total_weight(m_zero_weight)
    ,m_nonzero_count_lower(0), m_nonzero_count_upper(0)
    ,m_total_weight_lower(m_zero_weight), m_total_weight_upper(m_zero_weight)
#if DYNDIST_IFASSERT
    ,m_update_totals_needs_completion(false)
    ,m_batch_update_in_progress(false)
#endif
    ,m_levels()
    ,m_level_total(zero_level), m_level_tiny(zero_level)
    ,m_force_meta(0), m_meta(NULL)
#if DYNDIST_IFSTATS
    ,m_draws(0)
    ,m_draws_level_meta(0)
    ,m_draws_level_scan(0), m_draws_level_scan_iterations(0)
    ,m_draws_single(0)
    ,m_draws_acceptor(0), m_draws_acceptor_iterations(0)
    ,m_draws_scan(0), m_draws_scan_iterations(0)
    ,m_draws_purify(0)
#endif
{}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
discrete_distribution<S,W,D,C>::~discrete_distribution()
{
#if DYNDIST_IFVERIFY
    verify_consistency();
#endif
    if (m_meta)
        delete m_meta;
}

template<typename S,typename W,typename D,typename C>
template<typename Engine>
DYNDIST_NOINLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::operator()
(Engine& engine)
{
    DYNDIST_ASSERT(!m_batch_update_in_progress);
    
    if (m_nonzero_count == 0)
        throw std::out_of_range("distribution has total weight zero");
    DYNDIST_ASSERT(m_total_weight > m_zero_weight);
    
    typedef uniform_distribution<weight_type, Engine> uniform_t;
    typedef typename uniform_t::result_type x_t;
    x_t x = uniform_t::sample(m_total_weight, engine);
    return operator()(engine, x);
}

template<typename S,typename W,typename D,typename C>
template<typename Engine>
DYNDIST_NOINLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::operator()
(Engine& engine,
 typename uniform_distribution<weight_type, Engine>::result_type& x)
{
    DYNDIST_ASSERT(!m_batch_update_in_progress);
    DYNDIST_STATS(++m_draws);
    
    if (m_nonzero_count == 0)
        throw std::out_of_range("distribution has total weight zero");
    DYNDIST_ASSERT(m_total_weight > m_zero_weight);

    levels_value_type* level_ptr;
    if ((m_force_meta > 0) ||
         ((m_force_meta == 0) && (m_levels_nonzero.size() > meta_threshold)))
    {
        /* Use meta distribution */
        level_ptr = &draw_level_meta(engine, x);
    }
    else
    {
        /* Don't use meta distribution */
        level_ptr = &draw_level_scan(x);
    }
    
    levels_value_type& level = *level_ptr;
    const std::size_t level_size = level.second.events.size();
    if (level_size == 1) {
        DYNDIST_STATS(++m_draws_single);
        DYNDIST_ASSERT(x < level.second.total_weight);
        return pointer(level, level.second.events.front());
    }
    else if (!level.second.tiny && !level.second.pure) {
        return draw_purify(level, x);
    }
    else if (level.second.pure && (level_size > scan_threshold)) {
        return draw_acceptor(level, engine, x);
    }
    else {
        return draw_scan(level, x);
    }
}

template<typename S,typename W,typename D,typename C>
template<typename X>
DYNDIST_NOINLINE DYNDIST_FLATTEN
typename discrete_distribution<S,W,D,C>::levels_value_type&
discrete_distribution<S,W,D,C>::draw_level_scan
(X& x)
{
    DYNDIST_STATS(++m_draws_level_scan);
    
    DYNDIST_ASSERT(x < m_total_weight);
    
    /* Scan levels, looking for first cummulative weight to exceed <random> */
    typedef typename levels_collection_type::iterator iterator;
    iterator i = m_levels_nonzero.begin();
    const iterator end = m_levels_nonzero.end();
    DYNDIST_ASSERT(i != end);
    
    /* Test first element */
    DYNDIST_STATS(++m_draws_level_scan_iterations);
    if (x < (*i)->second.total_weight)
        return **i;
    
    /* Test second...last elements and do a bubble-sort pass */
    weight_type cummulative = (*i)->second.total_weight;
    for(iterator p=i++; i != end; ++i, ++p)
    {
        levels_value_type& level = **i;
        level_data_type& ld = level.second;
        
        /* Exchange entry with previous ones if previous weight was smaller.
         * Beware that i might afterwards point to a level that was already
         * processed! The code below thus uses "level" and "ld", and doesn't
         * dereference i again!
         */
        DYNDIST_ASSERT(p != end);
        if ((*p)->second.total_weight < ld.total_weight) {
            using std::swap;
            swap((*p)->second.link_nonzero, ld.link_nonzero);
            swap(*p, *i);
            DYNDIST_ASSERT(i - m_levels_nonzero.begin() ==
                         (std::ptrdiff_t)(*i)->second.link_nonzero);
            DYNDIST_ASSERT(p - m_levels_nonzero.begin() ==
                         (std::ptrdiff_t)(*p)->second.link_nonzero);
        }
        
        /* Test if current entry makes the cummulative weight exceed the sample */
        const weight_type cummulative_prev = cummulative;
        cummulative += ld.total_weight;
        DYNDIST_STATS(++m_draws_level_scan_iterations);
        if (x < cummulative) {
            x -= cummulative_prev;
            return level;
        }
    }
    
    throw std::logic_error("internal error, distribution is inconsistent");
}

template<typename S,typename W,typename D,typename C>
template<typename Engine>
DYNDIST_NOINLINE DYNDIST_FLATTEN
typename discrete_distribution<S,W,D,C>::levels_value_type&
discrete_distribution<S,W,D,C>::draw_level_meta
(Engine& engine,
 typename uniform_distribution<weight_type, Engine>::result_type& x)
{
    DYNDIST_STATS(++m_draws_level_meta);

    /*  Make sure the meta distribution is up-to-date */
    update_meta();
    
    /* Draw level from meta distribution */
    const meta_pointer meta_event = m_meta->operator()(engine, x);
    return *reinterpret_cast<levels_value_type*>(meta_event->data);
}

template<typename S,typename W,typename D,typename C>
template<typename Engine>
DYNDIST_NOINLINE DYNDIST_FLATTEN
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::draw_acceptor
(levels_value_type& level, Engine& engine,
 typename uniform_distribution<weight_type, Engine>::result_type& x)
{
    DYNDIST_STATS(++m_draws_acceptor);
    
    typedef pow2_expression<weight_type> pow2_expr_t;
    typedef std::pair<size_t, pow2_expr_t> bounds_t;
    typedef uniform_distribution<bounds_t, Engine> unif_sample_t;
    typedef typename unif_sample_t::result_type sample_t;

    const level_type lk = level.first;
    level_data_type& ld = level.second;
    DYNDIST_ASSERT(ld.pure);

    while (true) {
        DYNDIST_STATS(++m_draws_acceptor_iterations);
        
        /* Pick uniform sample from {0,...,N-1} x [0,W), where N is the number
         * of events in the level and W is the upper bound on the weight.
         */
        const bounds_t upper = bounds_t((size_t)ld.events.size(),
                                        pow2_expr_t(lk, ld.weight_upper));
        const sample_t sample = unif_sample_t::sample(upper, engine);
        const size_t sample_index = (size_t)sample.first;
        
        /* Use first component of the sample to pick an event */
        DYNDIST_ASSERT(sample_index < ld.events.size());
        value_type& event_value = ld.events[(std::size_t)sample_index];
        DYNDIST_ASSERT(event_value.weight <= ld.weight_upper);
        DYNDIST_ASSERT(ld.weight_lower < event_value.weight);
        
        /* Use to second component to accept or reject the event.
         * Since the weight is surely greather than weight_lower, a
         * r < weight_lower will accept any event. This is used to avoid
         * accessing the event at all for such r.
         */
        DYNDIST_ASSERT(sample.second < ld.weight_upper);
        if ((sample.second < ld.weight_lower) ||
            (sample.second < event_value.weight))
        {
            x = sample.second;
            return pointer(level, event_value);
        }
    }
}

template<typename S,typename W,typename D,typename C>
template<typename X>
DYNDIST_NOINLINE DYNDIST_FLATTEN
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::draw_scan
(levels_value_type& level, X& x)
{
    DYNDIST_STATS(++m_draws_scan);

    level_data_type& ld = level.second;
    DYNDIST_ASSERT(x < ld.total_weight);
    
    /* Scan events, looking for first cummulative weight to exceed <random> */
    weight_type cummulative = m_zero_weight;
    typedef typename level_events_type::iterator iterator_type;
    const iterator_type end = ld.events.end();
    for(iterator_type e = ld.events.begin(); e != end; ++e)
    {
        DYNDIST_STATS(++m_draws_scan_iterations);

        value_type& event_value = *e;
        DYNDIST_ASSERT(event_value.weight <= ld.weight_upper);

        const weight_type cummulative_prev = cummulative;
        cummulative += event_value.weight;
        if (x < cummulative) {
            x -= cummulative_prev;
            return pointer(level, event_value);
        }
    }
    
    throw std::logic_error("internal error, distribution is inconsistent");
}

template<typename S,typename W,typename D,typename C>
template<typename X>
DYNDIST_NOINLINE DYNDIST_FLATTEN
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::draw_purify
(levels_value_type& level, X& x)
{
    DYNDIST_STATS(++m_draws_purify);

    level_data_type& ld = level.second;
    DYNDIST_ASSERT(x < ld.total_weight);
    DYNDIST_ASSERT(!ld.pure);
    DYNDIST_ASSERT(!ld.tiny);
    
    /* Scan events, and
     *   (a) Look for first cummulative weight to exceed <random>
     *   (b) Move events to their correct level
     */
    
    const weight_type total_weight_before = ld.total_weight;
    weight_type cummulative = m_zero_weight;
    
    DYNDIST_ASSERT(!ld.events.empty());
    typedef typename level_events_type::iterator iterator_type;
    const iterator_type first = ld.events.begin();
    iterator_type last = ld.events.end() - 1;
    pointer result = pointer();
    for(iterator_type event_it = ld.events.end(); event_it != first; )
    {
        /* Dereference and *decrement* iterator */
        level_events_value_type& event_value_holder = *(--event_it);
        value_type& event_value = event_value_holder;
        DYNDIST_ASSERT(event_value.weight <= ld.weight_upper);
        DYNDIST_ASSERT(last == ld.events.end() - 1);
        
        /* (a) */
        if (!result) {
            /* No result yet */
            
            const weight_type cummulative_prev = cummulative;
            cummulative += event_value.weight;
            if (x < cummulative) {
                x -= cummulative_prev;
                result = pointer(level, event_value);
            }
            
            /* Note that the event that result now points to may be moved
             * before we return. The purification code below takes care of
             * updating the result if that happens.
             */
        }
#if DYNDIST_IFASSERT
        else {
            /* If asserts are enabled, continue to track the total weight */
            cummulative += event_value.weight;
        }
        DYNDIST_ASSERT(cummulative <= total_weight_before);
#endif
        
        /* (b) */
        if (event_value.weight <= ld.weight_lower) {
            /* Move event to its appropriate level */
            
            /* Copy event to the appropriate layer */
            const pointer event = pointer(level, event_value);
            const level_type l = weight_level(event_value.weight);
            DYNDIST_ASSERT(l < level.first);
            const pointer event_new = place(event_value, *m_levels.locate(l));
            
            /* Update result & notify */
            if (result == event)
                result = event_new;
            m_event_moved(event, event_new);
            
            /* Update level's total weight */
            ld.total_weight -= event_value.weight;
            
            /* Remove event by replacing it with the level's last entry */
            if (event_it != last) {
                value_type& last_value = *last;
                event_value_holder = last_value;
                
                /* Update result & notify */
                const pointer last_event = pointer(level, last_value);
                m_event_moved(last_event, event);
                if (result == last_event)
                    result = event;
            }
            
            /* Remove last entry */
            if (last != first) {
                last -= 1;
                ld.events.pop_back();
                DYNDIST_ASSERT(!ld.events.empty());
            }
            else {
                ld.events.pop_back();
                DYNDIST_ASSERT(ld.events.empty());
                break;
            }
        }
    }
    DYNDIST_ASSERT(cummulative == total_weight_before);
    
    /* Update level */
    ld.pure = true;
    if (ld.total_weight != total_weight_before) {
        levels_updated_add(level);
        if (ld.total_weight == m_zero_weight)
            levels_nonzero_remove(level);
    }
    
    if (!result)
        throw std::logic_error("failed to draw sample, distribution is inconsistent");
    
    return result;
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
void
discrete_distribution<S,W,D,C>::clear()
{
    DYNDIST_ASSERT(!m_batch_update_in_progress);
    m_zero_weight = weight_type();
    m_total_count = 0;
    m_nonzero_count = 0;
    m_total_weight = m_zero_weight;
    m_nonzero_count_lower = 0;
    m_nonzero_count_upper = 0;
    m_total_weight_lower = m_zero_weight;
    m_total_weight_upper = m_zero_weight;
#if DYNDIST_IFASSERT
    m_update_totals_needs_completion = false;
    m_batch_update_in_progress = false;
#endif
    m_levels.clear();
    m_level_zero = levels_value_type();
    m_levels_nonzero.clear();
    m_levels_updated.clear();
    m_level_total = zero_level;
    m_level_tiny = zero_level;
    if (m_meta) {
        delete m_meta;
        m_meta = NULL;
    }
#if DYNDIST_IFSTATS
    m_draws = 0;
    m_draws_level_meta = 0;
    m_draws_level_scan = 0;
    m_draws_level_scan_iterations = 0;
    m_draws_single = 0;
    m_draws_acceptor = 0;
    m_draws_acceptor_iterations = 0;
    m_draws_scan = 0;
    m_draws_scan_iterations = 0;
    m_draws_purify = 0;
#endif

#if DYNDIST_IFVERIFY
    verify_consistency();
#endif
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::insert
(const value_type& event)
{
    DYNDIST_ASSERT(!m_batch_update_in_progress);

    pointer result;
    
    /* Upon seeing the first user-supplied weight instance, make sure all the
     * weight instances we default constructed are replaced by instances that
     * are derived from the user-supplied one. This is necessary to support
     * flavoured weights.
     */
    if (m_total_count == 0)
        initialize_weight_type(event.weight);
    
    /* Increment total count */
    if (m_total_count >= max_size)
        throw std::length_error("maximal size exceeded");
    m_total_count += 1;
    
    /* Insert event */
    if (event.weight == m_zero_weight) {
        /* Place on zero level */
        result = place_zero(event);
    }
    else {
        /* Update totals */
        const bool utc = update_totals(m_nonzero_count + 1,
                                       m_total_weight + event.weight);
        
        /* Place on a non-zero level */
        result = place(event, *m_levels.locate(weight_level(event.weight)));
        
        /* Cleanup */
        if (utc)
            update_totals_complete();
    }
    
#if DYNDIST_IFVERIFY
    verify_consistency();
#endif
    
    return result;
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::batch_insert
(const value_type& event)
{
    /* Total weight must already have been updated here! */

    if (event.weight == m_zero_weight)
        return place_zero(event);
    else
        return place(event, *m_levels.locate(weight_level(event.weight)));
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::set
(const pointer event, const weight_type& weight)
{
    DYNDIST_ASSERT(!m_batch_update_in_progress);

    value_type& event_value = *event;
    const bool weight_zero = (weight == m_zero_weight);
    const bool weight_was_zero = (event_value.weight == m_zero_weight);
    
    /* Update totals */
    const bool utc = update_totals(m_nonzero_count
                                   + (weight_zero ? 0 : 1)
                                   - (weight_was_zero ? 0 : 1),
                                   m_total_weight + weight - event_value.weight);
    
    /* Update event's weight */
    const pointer result = update_weight(event, event_value,
                                         weight, weight_zero, weight_was_zero);
    
    /* Cleanup */
    if (utc)
        update_totals_complete();
    
#if DYNDIST_IFVERIFY
    verify_consistency();
#endif
    
    return result;
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::batch_set
(const pointer event, const weight_type& weight)
{
    /* Update existing element */
    value_type& event_value = *event;
    const bool weight_zero = (weight == m_zero_weight);
    const bool weight_was_zero = (event_value.weight == m_zero_weight);
    
    /* Total weight must already have been updated here! */
    
    return update_weight(event, event_value,
                         weight, weight_zero, weight_was_zero);
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
void
discrete_distribution<S,W,D,C>::erase
(const pointer event)
{
    DYNDIST_ASSERT(!m_batch_update_in_progress);

    const value_type& event_value = *event;

    /* Update totals */
    m_total_count -= 1;
    DYNDIST_ASSERT(m_total_count >= 0);
    
    if (event_value.weight == m_zero_weight) {
        /* Remove event from the zero level */
        unplace_zero(event);
    }
    else {
        /* Update totals */
        const bool utc = update_totals(m_nonzero_count - 1,
                                       m_total_weight - event_value.weight);
        
        /* Remove event from its non-zero level */
        unplace(event);
        
        /* Cleanup */
        if (utc)
            update_totals_complete();
    }
    
#if DYNDIST_IFVERIFY
    verify_consistency();
#endif
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
void
discrete_distribution<S,W,D,C>::batch_erase
(const pointer event)
{
    const value_type& event_value = *event;
    
    /* Remove event */
    if (event_value.weight == m_zero_weight)
        unplace_zero(event);
    else
        unplace(event);
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
discrete_distribution<S,W,D,C>::batch_update::batch_update
(discrete_distribution& distribution,
 const size_t total_count_final, const size_t nonzero_count_final,
 const weight_type& weight_final)
    :m_distribution(distribution)
#if DYNDIST_IFASSERT
    ,m_total_count(m_distribution.m_total_count)
    ,m_nonzero_count(m_distribution.m_nonzero_count)
    ,m_total_weight(m_distribution.m_total_weight)
#endif
{
    /* Upon seeing the first user-supplied weight instance, make sure all the
     * weight instances we default constructed are replaced by instances that
     * are derived from the user-supplied one. This is necessary to support
     * flavoured weights.
     */
    if (m_distribution.m_total_count == 0)
        m_distribution.initialize_weight_type(weight_final);
    
    /* Update totals */
    m_distribution.m_total_count = total_count_final;
    m_utc = m_distribution.update_totals(nonzero_count_final, weight_final);
    
#if DYNDIST_IFASSERT
    DYNDIST_ASSERT(!m_distribution.m_batch_update_in_progress);
    m_distribution.m_batch_update_in_progress = true;
#endif
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
discrete_distribution<S,W,D,C>::batch_update::~batch_update()
{
    /* Cleanup if necessary */
    if (m_utc)
        m_distribution.update_totals_complete();
    
#if DYNDIST_IFASSERT
    DYNDIST_ASSERT(m_total_count == m_distribution.m_total_count);
    DYNDIST_ASSERT(m_nonzero_count == m_distribution.m_nonzero_count);
    DYNDIST_ASSERT(m_total_weight == m_distribution.m_total_weight);
    DYNDIST_ASSERT(m_distribution.m_batch_update_in_progress);
    m_distribution.m_batch_update_in_progress = false;
#endif

#if DYNDIST_IFVERIFY
    m_distribution.verify_consistency();
#endif
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::batch_update::insert
(const value_type& event)
{
#if DYNDIST_IFASSERT
    m_total_count += 1;
    m_nonzero_count += (event.weight != m_distribution.m_zero_weight) ? 1 : 0;
    m_total_weight += event.weight;
#endif

    return m_distribution.batch_insert(event);
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::batch_update::set
(const pointer event, const weight_type& weight)
{
#if DYNDIST_IFASSERT
    m_nonzero_count += (((weight != m_distribution.m_zero_weight) ? 1 : 0) -
                        ((event->weight != m_distribution.m_zero_weight) ? 1 : 0));
    m_total_weight += weight;
    m_total_weight -= event->weight;
#endif

    return m_distribution.batch_set(event, weight);
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
void
discrete_distribution<S,W,D,C>::batch_update::erase
(const pointer event)
{
#if DYNDIST_IFASSERT
    m_total_count -= 1;
    m_nonzero_count -= (event->weight != m_distribution.m_zero_weight) ? 1 : 0;
    m_total_weight -= event->weight;
#endif
    
    return m_distribution.batch_erase(event);
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::update_weight
(const pointer event, value_type& event_value,
 const weight_type& weight,
 const bool weight_zero, const bool weight_was_zero)
{
    /* Total weight must already have been updated here! */
    
    if (!weight_was_zero && !weight_zero)
        return update_weight_nonzero_nonzero(event, weight);
    else if (!weight_was_zero)
        return update_weight_nonzero_zero(event);
    else if (!weight_zero)
        return update_weight_zero_nonzero(event, weight);
    else
        return event;
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::update_weight_nonzero_nonzero
(const pointer event, const weight_type& weight)
{
    value_type& event_value = *event;
    DYNDIST_ASSERT(weight > m_zero_weight);
    DYNDIST_ASSERT(event_value.weight > m_zero_weight);

    /* Total weight must already have been updated here! */
    
    /* Move to correct level */
    levels_value_type& level = *event.level;
    level_data_type& ld = level.second;
    if ((ld.weight_lower < weight) && (weight <= ld.weight_upper))
    {
        /* Level stayed the same */
        
        /* Update totals */
        DYNDIST_ASSERT(log2ceil(weight) == level.first);
        ld.total_weight -= event_value.weight;
        ld.total_weight += weight;
        DYNDIST_ASSERT(ld.total_weight > m_zero_weight);
        
        /* Update event */
        const_cast<weight_type&>(event_value.weight) = weight;

        /* Mark level as updated */
        levels_updated_add(level);
        
        return event;
    }
    else if (ld.tiny && (weight <= ld.weight_upper))
    {
        /* Level stayed the same, but only because its the tiny level */
        
        /* Update totals */
        DYNDIST_ASSERT(log2ceil(weight) < level.first);
        ld.total_weight -= event_value.weight;
        ld.total_weight += weight;
        ld.pure = false;
        DYNDIST_ASSERT(ld.total_weight > m_zero_weight);

        /* Update event */
        const_cast<weight_type&>(event_value.weight) = weight;

        /* Mark level as updated */
        levels_updated_add(level);
        
        return event;
    }
    else {
        /* Level changed */
        
        /* Create new event */
        const pointer event_new =
        place(value_type(*event, weight),
              *m_levels.locate(weight_level(weight)));
        
        /* Notify */
        m_event_moved(event, event_new);
        
        /* Remove old event */
        unplace(event);
        
        return event_new;
    }
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::update_weight_zero_nonzero
(const pointer event, const weight_type& weight)
{
    const value_type& event_value = *event;
    DYNDIST_ASSERT(event_value.weight == m_zero_weight);
    
    /* Total weight must already have been updated here! */
    
    /* Create new event */
    const pointer event_new = place(value_type(event_value, weight),
                                    *m_levels.locate(weight_level(weight)));
    
    /* Notify */
    m_event_moved(event, event_new);
    
    /* Remove old event */
    unplace_zero(event);
    
    return event_new;
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::update_weight_nonzero_zero
(const pointer event)
{
    const value_type& event_value = *event;
    DYNDIST_ASSERT(event->weight > m_zero_weight);
    
    /* Total weight may, but doesn't have to have beeen already updated here  */
    
    /* Create new event */
    const pointer event_new = place_zero(value_type(event_value, m_zero_weight));
    
    /* Notify */
    m_event_moved(event, event_new);
    
    /* Remove old event */
    unplace(event);
    
    return event_new;
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::place
(const value_type& event_value, levels_value_type& level)
{
    DYNDIST_ASSERT(event_value.weight > m_zero_weight);
    
    /* Level will now have non-zero weight */
    level_data_type& ld = level.second;
    if (ld.events.empty()) {
        DYNDIST_ASSERT(ld.total_weight == m_zero_weight);
        DYNDIST_ASSERT(ld.link_nonzero == -std::size_t(1));
        m_levels_nonzero.push_back(&level);
        ld.link_nonzero = m_levels_nonzero.size() - 1;
    }
    
    /* Validate that level is appropriate for weight */
    DYNDIST_ASSERT(event_value.weight <= ld.weight_upper);
    DYNDIST_ASSERT(log2ceil(event_value.weight) <= level.first);
    DYNDIST_ASSERT((ld.weight_lower < event_value.weight) || ld.tiny);
    DYNDIST_ASSERT((log2ceil(event_value.weight) == level.first) || ld.tiny);

    /* Add to level */
    ld.total_weight += event_value.weight;
    ld.events.push_back(event_value);
    const pointer event = pointer(level, ld.events.back());
    ld.pure = ld.pure && (ld.weight_lower < event_value.weight);
    
    /* Mark level as updated */
    levels_updated_add(level);
    
    return event;
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
typename discrete_distribution<S,W,D,C>::pointer
discrete_distribution<S,W,D,C>::place_zero
(const value_type& event_value)
{
    DYNDIST_ASSERT(event_value.weight == m_zero_weight);
    
    /* Add to zero level */
    DYNDIST_ASSERT(event_value.weight == m_zero_weight);
    m_level_zero.second.events.push_back(event_value);
    const pointer event = pointer(m_level_zero,
                                  m_level_zero.second.events.back());
    
    return event;
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
void
discrete_distribution<S,W,D,C>::unplace
(const pointer event)
{
    level_events_value_type& event_value_holder = event.as_value_holder();
    value_type& event_value = event_value_holder;
    DYNDIST_ASSERT(event_value.weight > m_zero_weight);
    
    /* Get event's current level */
    levels_value_type& level = *event.level;
    level_data_type& ld = level.second;
    DYNDIST_ASSERT(!ld.events.empty());
    DYNDIST_ASSERT(event_value.weight <= ld.weight_upper);
    DYNDIST_ASSERT(log2ceil(event_value.weight) <= level.first);
    
    /* Update level's total weight */
    ld.total_weight -= event_value.weight;
    levels_updated_add(level);
    if (ld.total_weight == m_zero_weight)
        levels_nonzero_remove(level);
    
    /* Replace event with last entry
     * Note that we *musn't* swap & notify if the element is the last element!
     * The caller of unplace() already notified about a new pointer for that
     * element, so notifying again with the old pointer would cause confusion
     */
    value_type& last = ld.events.back();
    if (&event_value != &last) {
        event_value_holder = last;
    
        /* Notify */
        m_event_moved(pointer(level, last), event);
    }
    
    /* Remove last entry */
    ld.events.pop_back();
    DYNDIST_ASSERT(ld.events.empty() == (ld.total_weight == m_zero_weight));
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
void
discrete_distribution<S,W,D,C>::unplace_zero
(const pointer event)
{
    level_events_value_type& event_value_holder = event.as_value_holder();
    value_type& event_value = event_value_holder;
    DYNDIST_ASSERT(event_value.weight == m_zero_weight);
    DYNDIST_ASSERT(event.level == &m_level_zero);
    
    /* Replace event with last entry
     * Note that we *musn't* swap & notify if the element is the last element!
     * The caller of unplace() already notified about a new pointer for that
     * element, so notifying again with the old pointer would cause confusion
     */
    value_type& last = m_level_zero.second.events.back();
    if (&event_value != &last) {
        event_value_holder = last;
        
        /* Notify */
        m_event_moved(pointer(m_level_zero, last), event);
    }
    
    /* Remove last entry */
    m_level_zero.second.events.pop_back();
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
typename discrete_distribution<S,W,D,C>::level_type
discrete_distribution<S,W,D,C>::weight_level
(const weight_type& weight)
{
    return std::max(log2ceil(weight), m_level_tiny);
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
void
discrete_distribution<S,W,D,C>::initialize_weight_type
(const weight_type& weight)
{
    /* Initialize weight-related things upon seeing the first actual weight */
    DYNDIST_ASSERT(m_levels.empty());
    DYNDIST_ASSERT(m_total_count == 0);
    DYNDIST_ASSERT(m_nonzero_count == 0);
    m_zero_weight = weight;
    m_zero_weight = 0;
    m_total_weight = m_zero_weight;
    m_total_weight_lower = m_zero_weight;
    m_total_weight_upper = m_zero_weight;
    m_level_zero.second.total_weight = m_zero_weight;
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
bool
discrete_distribution<S,W,D,C>::update_totals
(size_t nonzero_count, const weight_type& total_weight)
{
    m_total_weight = total_weight;
    m_nonzero_count = nonzero_count;
    
    /* Fast path if bounds stay within slack we accomodated for */
    const bool weight_within_slack = ((m_total_weight_lower < m_total_weight) &&
                                      (m_total_weight <= m_total_weight_upper));
    const bool count_within_slack = ((m_nonzero_count_lower < m_nonzero_count) &&
                                     (m_nonzero_count <= m_nonzero_count_upper));
    if (weight_within_slack && count_within_slack)
    {
        DYNDIST_ASSERT(!m_levels.empty());
        DYNDIST_ASSERT(nonzero_count > 0);
        
        DYNDIST_ASSERT(m_level_tiny <= (log2ceil(m_total_weight) - 1 -
                                      2*log2ceil(m_nonzero_count)));
        DYNDIST_ASSERT(m_level_total >= log2ceil(m_total_weight));
        
        /* If there are no uncomplete totals updates, match must be precise */
#if DYNDIST_IFASSERT
        if (!m_update_totals_needs_completion) {
            DYNDIST_ASSERT(m_level_tiny == m_levels.lower_bound());
            DYNDIST_ASSERT(m_level_total == m_levels.upper_bound());
        }
        else {
            DYNDIST_ASSERT(m_level_tiny >= m_levels.lower_bound());
            DYNDIST_ASSERT(m_level_total <= m_levels.upper_bound());
        }
#endif
        return false;
    }
    
    /* Slow path */
    update_totals_begin(weight_within_slack, count_within_slack);
    return true;
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
void
discrete_distribution<S,W,D,C>::update_totals_begin
(const bool weight_within_slack, const bool count_within_slack)
{
    /* Update "tiny" and "total" thresholds.
     *
     * The maximum weight of a single event is obviously at most total_weight.
     * It thus suffices for levels up to <total>, where
     *   total >= log2ceil(total_weight)
     * to exist.
     *
     * On the other end of the scale, the goal is to put all weights which are
     * no larger than total_weight / count^2 into the tiny level -- that ensures
     * that the cummulative probability of the tiny level is at most 1/count^2,
     * meaning that handling it with a linear-time algorithmn keeps the overall
     * expected selection time O(1). This is ensured by putting all events with
     *   log2ceil(weight) <= tiny
     * into the tiny level, where
     *   tiny <= log2(total_weight / nonzero_count^2).
     *
     * We avoid repeated recomputations of tiny and total by allowing for some
     * slack. Whenever we *do* have to recompute, we compute a value of tiny
     * and total which works for all values of total_weight within
     *   (total_weight_lower, total_weight_upper],
     * and for all values of nonzero_count within
     *   (nonzero_count_lower, nonzero_count_upper],
     * i.e. we set
     *   tiny = log2(total_weight_lower) - log2(nonzero_count_upper).
     *
     * The interval bounds are always integrals powers of two, and chosen such
     * that repeated recomputations are still avoided even if one or both of
     * the values oscillated around some integral power of two. Assuming that
     * w = log2ceil(total_weight), and c = log2ceil(nonzero_count), we pick
     *   (total_weight_lower, total_weight_upper] = (2^(w-2), 2^w],
     *   (nonzero_count_lower , nonzero_count_upper ] = (2^(c-2), 2^c].
     * This ensures that whenver total_weight and nonzero_count lie within
     * those intervals, then
     *   w-2 <= log2(total_weight), and
     *   log2(nonzero_count) <= c
     */

    /* Clear all non-zero levels if the last non-zero weight is removed */
    if (m_nonzero_count == 0) {
        m_total_weight_lower = m_zero_weight;
        m_total_weight_upper = m_zero_weight;
        m_nonzero_count_lower = 0;
        m_nonzero_count_upper = 0;
        m_level_tiny = zero_level;
        m_level_total = zero_level;
#if DYNDIST_IFASSERT
        m_update_totals_needs_completion = true;
#endif
        return;
    }

    /* Clear "tiny" flag of current tiny level */
    if (!m_levels.empty())
        m_levels[m_level_tiny].tiny = false;
    
    /* Handle changed total_weight */
    const level_type SL = 2;
    if (!weight_within_slack) {
        const level_type total_weight_log2ceil = log2ceil(m_total_weight);
        m_total_weight_lower = pow2(total_weight_log2ceil - SL, m_zero_weight);
        m_total_weight_upper = pow2(total_weight_log2ceil    , m_zero_weight);
        m_level_tiny += (total_weight_log2ceil - m_level_total);
        m_level_total = total_weight_log2ceil;
    }
    DYNDIST_ASSERT((m_total_weight_lower < m_total_weight) &&
                 (m_total_weight <= m_total_weight_upper));
    
    /* Handle changed nonzero_count */
    if (!count_within_slack) {
        const level_type nonzero_count_log2ceil = log2ceil(m_nonzero_count);
        m_nonzero_count_lower = pow2(nonzero_count_log2ceil - SL, size_t(0));
        m_nonzero_count_upper = pow2(nonzero_count_log2ceil    , size_t(0));
        m_level_tiny = m_level_total - SL - 2*nonzero_count_log2ceil;
    }
    DYNDIST_ASSERT((m_nonzero_count_lower < m_nonzero_count) &&
                 (m_nonzero_count <= m_nonzero_count_upper));

    /* Validate the end result */
    DYNDIST_ASSERT(m_level_tiny <= (log2ceil(m_total_weight) - 1 -
                                  2*log2ceil(m_nonzero_count)));
    DYNDIST_ASSERT(m_level_total >= log2ceil(m_total_weight));

    /* Expand level range */
    f_level_init level_init = f_level_init(*this);
    m_levels.enlarge_to(m_level_tiny, m_level_total, level_init);
    DYNDIST_ASSERT(!m_meta || (m_meta->size() == m_levels.size()));
    
    /* Mark new tiny level as "tiny" */
    m_levels[m_level_tiny].tiny = true;
    
    /* Levels outside [m_level_tiny, m_level_total] must be removed later  */
#if DYNDIST_IFASSERT
    m_update_totals_needs_completion = true;
#endif
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
void
discrete_distribution<S,W,D,C>::update_totals_complete()
{
    DYNDIST_ASSERT(m_update_totals_needs_completion);
    
    if (m_nonzero_count == 0) {
        /* Remove all non-zero levels and the meta-distribution */
        DYNDIST_ASSERT(m_total_weight == m_zero_weight);
        DYNDIST_ASSERT(m_total_count == m_level_zero.second.events.size());
#if DYNDIST_IFASSERT
        m_update_totals_needs_completion = false;
#endif
        m_levels.clear();
        m_levels_nonzero.clear();
        m_levels_updated.clear();
        m_level_tiny = zero_level;
        m_level_total = zero_level;
        if (m_meta) {
            delete m_meta;
            m_meta = NULL;
        }
        return;
    }

    /* Merge levels below the tiny level into the tiny level */
    levels_value_type& tiny = *m_levels.locate(m_level_tiny);
    for(typename levels_type::iterator i = m_levels.begin();
        i->first < m_level_tiny;
        m_levels.pop_front(), i = m_levels.begin())
    {
        DYNDIST_ASSERT(i != m_levels.end());

        /* Nothing to do for empty levels */
        levels_value_type& level = *i;
        level_data_type& ld = level.second;
        if (!ld.events.empty()) {
            /* Move events */
            const typename level_events_type::iterator end = ld.events.end();
            for(typename level_events_type::iterator e = ld.events.begin();
                e != end;
                ++e)
            {
                /* Event must belong to a smaller level */
                const value_type& event = *e;
                DYNDIST_ASSERT(event.weight <= tiny.second.weight_upper);
                
                /* Copy to the tiny layer */
                const pointer event_new = place(event, tiny);

                /* Notify */
                m_event_moved(pointer(level, *e), event_new);
            }
            
            /* Remove level from the non-zero level list */
            DYNDIST_ASSERT(ld.total_weight > m_zero_weight);
            ld.total_weight = m_zero_weight;
            levels_nonzero_remove(level);
        }
        else {
            DYNDIST_ASSERT(ld.total_weight == m_zero_weight);
            DYNDIST_ASSERT(ld.link_nonzero == -std::size_t(1));
        }
    
        /* Remove from the updated list */
        levels_updated_remove(level);

        /* If we have a meta distribution, remove level there too */
        if (m_meta)
            m_meta->erase(ld.link_meta);
    }
    DYNDIST_ASSERT(m_levels.front().second.tiny == true);
    
    /* Remove levels beyond level_total. Note that these aren't necessarily
     * empty -- some of these levels might have served as the tiny level
     * in the past, and so may contain events which weight_levels less than
     * the level's native weight level. Such events are moved to their correct
     * level.
     */
    for(typename levels_type::reverse_iterator i = m_levels.rbegin();
        i->first > m_level_total;
        m_levels.pop_back(), i = m_levels.rbegin())
    {
        DYNDIST_ASSERT(i != m_levels.rend());

        /* Nothing to do for empty levels */
        levels_value_type& level = *i;
        level_data_type& ld = level.second;
        if (!ld.events.empty()) {
            /* Move events */
            typename level_events_type::iterator end = ld.events.end();
            for(typename level_events_type::iterator e = ld.events.begin();
                e != end;
                ++e)
            {
                /* Events weight mustn't exceed total weight */
                const value_type& event = *e;
                DYNDIST_ASSERT(event.weight <= m_total_weight);
                
                /* Copy to the tiny layer */
                const pointer event_new =
                place(event, *m_levels.locate(weight_level(event.weight)));
                
                /* Notify */
                m_event_moved(pointer(level, *e), event_new);
            }
            
            /* Remove level from the non-zero level list */
            DYNDIST_ASSERT(ld.total_weight > m_zero_weight);
            ld.total_weight = m_zero_weight;
            levels_nonzero_remove(level);
        }
        else {
            DYNDIST_ASSERT(ld.total_weight == m_zero_weight);
            DYNDIST_ASSERT(ld.link_nonzero == -std::size_t(1));
        }
        
        /* Remove from the updated list */
        levels_updated_remove(level);
        
        /* If we have a meta distribution, remove level there too */
        if (m_meta)
            m_meta->erase(ld.link_meta);
    }

#if DYNDIST_IFASSERT
    DYNDIST_ASSERT(!m_meta || (m_meta->size() == m_levels.size()));
    m_update_totals_needs_completion = false;
#endif
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
void
discrete_distribution<S,W,D,C>::levels_updated_add
(levels_value_type& level)
{
    level_data_type& ld = level.second;

    if (ld.link_updated == -std::size_t(1)) {
        m_levels_updated.push_back(&level);
        ld.link_updated = m_levels_updated.size() - 1;
    }
    else {
        DYNDIST_ASSERT(m_levels_updated[ld.link_updated] == &level);
    }
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
void
discrete_distribution<S,W,D,C>::levels_updated_remove
(levels_value_type& level)
{
    level_data_type& ld = level.second;

    /* If the level is on the updated list, remove it */
    if (ld.link_updated != -std::size_t(1)) {
        /* Replace level's entry with the last entry */
        levels_value_type*& ud = m_levels_updated[ld.link_updated];
        DYNDIST_ASSERT(ud == &level);
        levels_value_type*& ud_last = m_levels_updated.back();
        ud = ud_last;
        m_levels_updated.pop_back();
        
        /* Update last entry's backlink */
        levels_value_type& level_last = *ud_last;
        level_last.second.link_updated = ld.link_updated;
        
        /* Reset level's nonzero-link */
        ld.link_updated = -std::size_t(1);
    }
}

template<typename S,typename W,typename D,typename C>
DYNDIST_INLINE
void
discrete_distribution<S,W,D,C>::levels_nonzero_remove
(levels_value_type& level)
{
    level_data_type& ld = level.second;
    
    /* Replace level's entry with the last entry */
    levels_value_type*& nz = m_levels_nonzero[ld.link_nonzero];
    DYNDIST_ASSERT(nz == &level);
    levels_value_type*& nz_last = m_levels_nonzero.back();
    nz = nz_last;
    m_levels_nonzero.pop_back();
    
    /* Update last entry's backlink */
    levels_value_type& level_last = *nz_last;
    level_last.second.link_nonzero = ld.link_nonzero;
    
    /* Reset level's nonzero-link */
    ld.link_nonzero = -std::size_t(1);
}

template<typename S,typename W,typename D,typename C>
DYNDIST_NOINLINE
void
discrete_distribution<S,W,D,C>::update_meta()
{
    DYNDIST_ASSERT(m_force_meta >= 0);
    
    if (!m_meta) {
        /* No meta distribution yet, create one */
        m_meta = new meta_type();

        /* Start batch update of meta distribution */
        typename meta_type::batch_update batch(*m_meta,
                                               (size_t)m_levels.size(),
                                               (size_t)m_levels_nonzero.size(),
                                               m_total_weight);

        /* And create events for all levels */
        DYNDIST_ASSERT(!m_levels.empty());
        const typename levels_type::iterator end = m_levels.end();
        for(typename levels_type::iterator i = m_levels.begin();
            i != end;
            ++i)
        {
            levels_value_type& level = *i;
            level_data_type& ld = level.second;

            ld.link_meta = batch.insert(meta_value_type(ld.total_weight, &level));
        }
    }
    else {
        /* Start batch update of meta distribution */
        typename meta_type::batch_update batch(*m_meta,
                                               (size_t)m_levels.size(),
                                               (size_t)m_levels_nonzero.size(),
                                               m_total_weight);
        
        /* Update meta distribution weights for modified levels */
        const typename levels_collection_type::const_iterator end = m_levels_updated.end();
        for(typename levels_collection_type::const_iterator i = m_levels_updated.begin();
            i != end;
            ++i)
        {
            levels_value_type& level = **i;
            level_data_type& ld = level.second;
            batch.set(ld.link_meta, ld.total_weight);
            ld.link_updated = -std::size_t(1);
        }
        m_levels_updated.clear();
    }
    
#if DYNDIST_IFVERIFY
    verify_meta_consistency();
#endif
}

#if DYNDIST_IFSTATS

template<typename S,typename W,typename D,typename C>
void
discrete_distribution<S,W,D,C>::print_statistics(std::ostream& dst,
                                            const std::size_t meta_depth)
const
{
    if (meta_depth == 0)
        dst << "*** Distribution Statistics ***\n";
    else
        dst << "*** Meta-" << meta_depth << " Distribution Statistics ***\n";
    dst << "| Elements with positive Weight           : ";
    dst << m_nonzero_count << '\n';
    dst << "| Levels                                  : ";
    dst << m_levels.size() << '\n';
    dst << "| Levels with positive Weight             : ";
    dst << m_levels_nonzero.size() << '\n';
    dst << "| Total Draws                             : ";
    dst << m_draws << '\n';
    dst << std::fixed << std::setprecision(8);
    if (m_draws > 0) {
        dst << "| | % which drew level from meta          : ";
        dst << double(m_draws_level_meta)/double(m_draws) << '\n';
        dst << "| | % which drew level using scan         : ";
        dst << double(m_draws_level_scan)/double(m_draws) << '\n';
        if (m_draws_level_scan > 0) {
            dst << "| | | average iterations required         : ";
            dst << double(m_draws_level_scan_iterations)/double(m_draws_level_scan) << '\n';
        }
        dst << "| | % which drew from singleton level     : ";
        dst << double(m_draws_single)/double(m_draws) << '\n';
        dst << "| | % which drew from level with acceptor : ";
        dst << double(m_draws_acceptor)/double(m_draws) << '\n';
        if (m_draws_acceptor > 0) {
            dst << "| | | average iterations required         : ";
            dst << double(m_draws_acceptor_iterations)/double(m_draws_acceptor) << '\n';
        }
        dst << "| | % which drew from level with scan     : ";
        dst << double(m_draws_scan)/double(m_draws) << '\n';
        if (m_draws_scan > 0) {
            dst << "| | | average iterations required         : ";
            dst << double(m_draws_scan_iterations)/double(m_draws_scan) << '\n';
        }
        dst << "| | % which drew from level and purified  : ";
        dst << double(m_draws_purify)/double(m_draws) << '\n';
    }
    if (m_meta)
        m_meta->print_statistics(dst, meta_depth + 1);
}

#endif /* DYNDIST_IFSTATS */

#if DYNDIST_IFVERIFY

template<typename S,typename W,typename D,typename C>
void
discrete_distribution<S,W,D,C>::verify_consistency() const
{
    /* Verify that all actions have been completed */
    DYNDIST_ASSERT(!m_batch_update_in_progress);
    DYNDIST_ASSERT(!m_update_totals_needs_completion);
    
    /* Handle empty domain */
    DYNDIST_ASSERT(m_zero_weight + m_zero_weight == m_zero_weight);
    if (m_levels.empty()) {
        DYNDIST_ASSERT(m_total_weight == m_zero_weight);
        DYNDIST_ASSERT(m_nonzero_count == 0);
        DYNDIST_ASSERT(m_levels.empty());
        DYNDIST_ASSERT(m_levels_updated.empty());
        DYNDIST_ASSERT(m_level_tiny == zero_level);
        DYNDIST_ASSERT(m_level_total == zero_level);
        DYNDIST_ASSERT(m_meta == NULL);
        return ;
    }
    else if ((m_total_weight == m_zero_weight) || (m_nonzero_count == 0)) {
        DYNDIST_ASSERT(m_total_weight == m_zero_weight);
        DYNDIST_ASSERT(m_nonzero_count == 0);
        DYNDIST_ASSERT(m_levels.empty());
        DYNDIST_ASSERT(m_levels_updated.empty());
        DYNDIST_ASSERT(m_level_tiny == zero_level);
        DYNDIST_ASSERT(m_level_total == zero_level);
        DYNDIST_ASSERT(m_meta == NULL);
        return;
    }
    DYNDIST_ASSERT(m_total_weight > m_zero_weight);
    DYNDIST_ASSERT(m_nonzero_count > 0);
    
    /* Validate order of levels */
    DYNDIST_ASSERT(zero_level < m_levels.lower_bound());
    DYNDIST_ASSERT(m_levels.lower_bound() == m_level_tiny);
    DYNDIST_ASSERT(m_level_tiny < m_level_total);
    DYNDIST_ASSERT(m_level_total == m_levels.upper_bound());
    
    /* Validate total counts and weight */
    weight_type total_weight = m_zero_weight;
    size_t nonzero_count = 0;
    
    /* Scan levels */
    for(typename levels_type::const_iterator i = m_levels.begin();
        i != m_levels.end();
        ++i)
    {
        /* Validate levels entry */
        const levels_value_type& level = *i;
        const level_type lk = i->first;
        DYNDIST_ASSERT(lk - m_levels.lower_bound() == i - m_levels.begin());
        const level_data_type& ld = i->second;
        DYNDIST_ASSERT(lk >= m_level_tiny - 1);
        DYNDIST_ASSERT(lk <= m_level_total + 1);
        
        /* First level should be marked "tiny" */
        DYNDIST_ASSERT(ld.tiny == (i == m_levels.begin()));
        
        /* Validate backlink into non-zero levels list */
        if (!ld.events.empty()) {
            DYNDIST_ASSERT(ld.total_weight > m_zero_weight);
            DYNDIST_ASSERT(ld.link_nonzero < m_levels_nonzero.size());
            DYNDIST_ASSERT(m_levels_nonzero[ld.link_nonzero] == &level);
        }
        else {
            DYNDIST_ASSERT(ld.total_weight == m_zero_weight);
            DYNDIST_ASSERT(ld.link_nonzero == -std::size_t(1));
        }

        /* Validate backlink into updated levels list */
        if (ld.link_updated != -std::size_t(1)) {
            DYNDIST_ASSERT(ld.link_updated < m_levels_updated.size());
            DYNDIST_ASSERT(m_levels_updated[ld.link_updated] == &level);
        }

        /* Validate weight range
         * Note that, since we don't check when creating levels whether weights
         * that belong to them are actually representable by weight_type, we
         * may have levels with a zero lower bound, or both a zero upper and
         * a zero lower bound
         */
        if (ld.weight_upper == m_zero_weight) {
            DYNDIST_ASSERT(ld.weight_lower == m_zero_weight);
            DYNDIST_ASSERT(ld.events.empty());
            DYNDIST_ASSERT(ld.total_weight == m_zero_weight);
        }
        else {
            DYNDIST_ASSERT(ld.weight_upper == pow2(lk, ld.total_weight));
            DYNDIST_ASSERT(log2ceil(ld.weight_upper) == lk);
            if (ld.weight_lower != m_zero_weight) {
                DYNDIST_ASSERT(ld.weight_lower == pow2(lk - 1, ld.total_weight));
                DYNDIST_ASSERT(log2ceil(ld.weight_lower) == lk - 1);
            }
        }

        /* Scan level's events */
        bool pure = true;
        weight_type level_weight = m_zero_weight;
        for(typename level_events_type::const_iterator e =
            ld.events.begin();
            e != ld.events.end();
            ++e)
        {
            /* Validate entry */
            const value_type& event = *e;
            DYNDIST_ASSERT(event.weight > m_zero_weight);
            
            /* Track count and total weight */
            nonzero_count += 1;
            level_weight += event.weight;

            /* Verify that the event belongs on this level */
            DYNDIST_ASSERT(event.weight <= ld.weight_upper);
            DYNDIST_ASSERT(log2ceil(event.weight) <= lk);
            DYNDIST_ASSERT((event.weight > ld.weight_lower) || !ld.pure);
            DYNDIST_ASSERT((log2ceil(event.weight) == lk) || !ld.pure);
            pure = pure && (ld.weight_lower < event.weight);
        }
        DYNDIST_ASSERT(ld.total_weight == level_weight);
        DYNDIST_ASSERT(pure || !ld.pure);
        
        total_weight += level_weight;
    }
    
    /* Scan zero level events */
    for(typename level_events_type::const_iterator e =
        m_level_zero.second.events.begin();
        e != m_level_zero.second.events.end();
        ++e)
    {
        /* Validate entry */
        const value_type& event = *e;
        DYNDIST_ASSERT(event.weight == m_zero_weight);
        total_weight += event.weight;
    }
    
    /* Validate total counts and weight */
    DYNDIST_ASSERT(m_total_weight == total_weight);
    DYNDIST_ASSERT(m_nonzero_count == nonzero_count);
    DYNDIST_ASSERT(m_total_count == m_nonzero_count + m_level_zero.second.events.size());
    
    /* Scan non-zero levels list */
    for(typename levels_collection_type::const_iterator i = m_levels_nonzero.begin();
        i != m_levels_nonzero.end();
        ++i)
    {
        /* Validate entry */
        const levels_value_type& level = **i;
        const std::ptrdiff_t lk = level.first;
        const level_data_type& ld = level.second;
        DYNDIST_ASSERT(m_levels.contains(lk));
        DYNDIST_ASSERT(&*m_levels.locate(lk) == *i);

        /* Validate that weight is indeed non-zero */
        DYNDIST_ASSERT(ld.total_weight > m_zero_weight);
    }

    /* Scan updated levels list */
    for(typename levels_collection_type::const_iterator i = m_levels_updated.begin();
        i != m_levels_updated.end();
        ++i)
    {
        /* Validate entry */
        const levels_value_type& level = **i;
        const std::ptrdiff_t lk = level.first;
        DYNDIST_ASSERT(m_levels.contains(lk));
        DYNDIST_ASSERT(&*m_levels.locate(lk) == *i);
    }
}

template<typename S,typename W,typename D,typename C>
void
discrete_distribution<S,W,D,C>::verify_meta_consistency()
{
    /* Validate existance of meta distribution */
    DYNDIST_ASSERT(m_meta);
    DYNDIST_ASSERT(m_force_meta >= 0);
    
    /* Validate meta distribution bounds */
    DYNDIST_ASSERT(m_meta->size() == m_levels.size());
    
    /* Validate meta distribution's internal consistency */
    m_meta->verify_consistency();

    /* Scan meta events */
    for(typename meta_type::iterator meta_i = m_meta->begin();
        meta_i != m_meta->end();
        ++meta_i)
    {
        const levels_value_type& level = *reinterpret_cast<levels_value_type*>(meta_i->data);
        DYNDIST_ASSERT(&level == &*m_levels.find(level.first));
        
        /* Validate crosslinks */
        DYNDIST_ASSERT(level.second.link_meta == (typename meta_type::pointer)meta_i);
        
        /* Indices and weights must agree */
        DYNDIST_ASSERT(meta_i->weight == level.second.total_weight);
    }

    /* Scan levels */
    for(typename levels_type::iterator level_i =  m_levels.begin();
        level_i != m_levels.end();
        ++level_i)
    {
        const levels_value_type& level = *level_i;
        
        /* Validate crosslinks */
        DYNDIST_ASSERT(level.second.link_meta->data == &level);
        typename meta_type::iterator meta_i;
        for(meta_i = m_meta->begin();
            (typename meta_type::pointer)meta_i != level.second.link_meta;
            ++meta_i)
            DYNDIST_ASSERT(meta_i != m_meta->end());
        
        /* Indices and weights must agree */
        DYNDIST_ASSERT(level.second.link_meta->weight == level.second.total_weight);
    }
}

#endif /* DYNDIST_IFVERIFY */
    
DYNDIST_NAMESPACE_END

#endif
