/*
 * dyndist/macros.h
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
#ifndef dyndist_macros_h
#define dyndist_macros_h

#ifndef __has_feature
#   define __has_feature(x) 0
#endif

#ifndef __has_attribute
#   define __has_attribute(x) 0
#endif

#ifndef _Pragma
#   define _Pragma(x)
#endif

#if __has_feature(cxx_inline_namespaces)
#   define DYNDIST_NAMESPACE_BEGIN namespace dyndist { inline namespace DYNDIST_FLAVOUR {
#   define DYNDIST_NAMESPACE_END } }
#   define DYNDIST_NAMESPACE dyndist
#else
#   define DYNDIST_NAMESPACE_BEGIN namespace dyndist { namespace DYNDIST_FLAVOUR {
#   define DYNDIST_NAMESPACE_END } using namespace DYNDIST_FLAVOUR; }
#   define DYNDIST_NAMESPACE dyndist::DYNDIST_FLAVOUR
#endif

#if __has_attribute(always_inline)
#   define DYNDIST_ATTR_INLINE __attribute__ ((always_inline))
#   ifndef DYNDIST_INLINE
#       define DYNDIST_INLINE inline DYNDIST_ATTR_INLINE
#   endif
#else
#   define DYNDIST_ATTR_INLINE
#   ifndef DYNDIST_INLINE
#       define DYNDIST_INLINE inline
#   endif
#endif

#if __has_attribute(flatten)
#   define DYNDIST_ATTR_FLATTEN __attribute__ ((flatten))
#else
#   define DYNDIST_ATTR_FLATTEN
#endif

#if __has_attribute(noinline)
#   define DYNDIST_ATTR_NOINLINE __attribute__ ((noinline))
#else
#   define DYNDIST_ATTR_NOINLINE
#endif

#ifndef DYNDIST_FLATTEN
#   define DYNDIST_FLATTEN DYNDIST_ATTR_FLATTEN
#endif
#ifndef DYNDIST_INLINE_FLATTEN
#   define DYNDIST_INLINE_FLATTEN DYNDIST_INLINE DYNDIST_FLATTEN
#endif
#ifndef DYNDIST_NOINLINE
#   define DYNDIST_NOINLINE DYNDIST_ATTR_NOINLINE
#endif

#if defined(__GNUC__)
#   define DYNDIST_NOWARN_PUSH \
    _Pragma("GCC diagnostic push") \
    _Pragma("GCC diagnostic ignored \"-Wunknown-pragmas\"")
#   define DYNDIST_NOWARN_SIGNS \
    _Pragma("GCC diagnostic ignored \"-Wsign-compare\"") \
    _Pragma("GCC diagnostic ignored \"-Wsign-conversion\"")
#   define DYNDIST_NOWARN_ALL \
    _Pragma("GCC diagnostic ignored \"-Wall\"") \
    _Pragma("GCC diagnostic ignored \"-Wshadow\"") \
    _Pragma("GCC diagnostic ignored \"-Wsign-compare\"") \
    _Pragma("GCC diagnostic ignored \"-Wsign-conversion\"") \
    _Pragma("GCC diagnostic ignored \"-Wconversion\"") \
    _Pragma("GCC diagnostic ignored \"-Wunused-parameter\"") \
    _Pragma("GCC diagnostic ignored \"-Wunused-variable\"") \
    _Pragma("GCC diagnostic ignored \"-Wdocumentation\"")
#   define DYNDIST_NOWARN_POP \
    _Pragma("GCC diagnostic pop")
#   define DYNDIST_EXPORT_BEGIN \
    _Pragma("GCC visibility push(default)")
#   define DYNDIST_EXPORT_END \
    _Pragma("GCC visibility pop")
#   define DYNDIST_NOEXPORT_BEGIN \
    _Pragma("GCC visibility push(hidden)")
#   define DYNDIST_NOEXPORT_END \
    _Pragma("GCC visibility pop")
#else
#   define DYNDIST_NOWARN_PUSH
#   define DYNDIST_NOWARN_SIGNS
#   define DYNDIST_NOWARN_ALL
#   define DYNDIST_NOWARN_POP
#   define DYNDIST_EXPORT_BEGIN
#   define DYNDIST_EXPORT_END
#   define DYNDIST_NOEXPORT_BEGIN
#   define DYNDIST_NOEXPORT_END
#endif

#define _DYNDIST_STRINGIFY(s) #s
#define DYNDIST_STRINGIFY(s) _DYNDIST_STRINGIFY(s)

#endif
