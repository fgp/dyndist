/*
 * dyndist/global.h
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

#ifndef dyndist_global_h
#define dyndist_global_h

#include "dyndist/macros.h"

#ifndef DYNDIST_ENABLE_VERIFICATION
#   define DYNDIST_ENABLE_VERIFICATION 0
#endif

#ifndef DYNDIST_ENABLE_ASSERTS
#   define DYNDIST_ENABLE_ASSERTS 0
#endif

#ifndef DYNDIST_ENABLE_INSTRUMENTATION
#   define DYNDIST_ENABLE_INSTRUMENTATION 0
#endif

#if DYNDIST_ENABLE_VERIFICATION
#   define DYNDIST_IFVERIFY 1
#   undef DYNDIST_ENABLE_ASSERTIONS
#   define DYNDIST_ENABLE_ASSERTIONS 1
#else
#   define DYNDIST_IFVERIFY 0
#endif

#if DYNDIST_ENABLE_ASSERTIONS
#   ifndef DYNDIST_ASSERT
#       include <cassert>
#       define DYNDIST_ASSERT assert
#   endif
#   define DYNDIST_IFASSERT 1
#   undef  DYNDIST_ENABLE_INSTRUMENTATION
#   define DYNDIST_ENABLE_INSTRUMENTATION 1
#else
#   undef  DYNDIST_ASSERT
#   define DYNDIST_ASSERT(x)
#   define DYNDIST_IFASSERT 0
#endif

#if DYNDIST_ENABLE_INSTRUMENTATION
#   define DYNDIST_STATS(x) x
#   define DYNDIST_IFSTATS 1
#else
#   define DYNDIST_STATS(x)
#   define DYNDIST_IFSTATS 0
#endif

#if DYNDIST_IFVERIFY
#   define DYNDIST_FLAVOUR vrf
#elif DYNDIST_IFASSERT
#   define DYNDIST_FLAVOUR chk
#elif DYNDIST_IFSTATS
#   define DYNDIST_FLAVOUR sts
#else
#   define DYNDIST_FLAVOUR dfl
#endif

#endif
