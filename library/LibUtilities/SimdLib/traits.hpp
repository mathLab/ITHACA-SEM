///////////////////////////////////////////////////////////////////////////////
//
// File: traits.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Vector type traits and tags.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITES_SIMDLIB_TRAITS_H
#define NEKTAR_LIB_LIBUTILITES_SIMDLIB_TRAITS_H

#include <type_traits>

namespace tinysimd
{

// Load tags
static constexpr struct is_aligned_t {} is_aligned{};
static constexpr struct is_not_aligned_t {} is_not_aligned{};
static constexpr struct is_not_reused_t {} is_not_reused{}; // streaming, skip cache

// Check load tags
template <class T>
struct is_streaming : std::integral_constant
    <
        bool,
        std::is_same<is_not_reused_t, typename std::remove_cv<T>::type>::value
    > {};

template <class T>
struct is_requiring_alignment : std::integral_constant
    <
        bool,
        std::is_same<is_aligned_t, typename std::remove_cv<T>::type>::value ||
        is_streaming<T>::value
    > {};

// Helper c++17 style
// template <class T>
// inline constexpr bool is_streaming_v = is_streaming<T>::value;
// template <class T>
// inline constexpr bool is_requiring_alignment_v = is_requiring_alignment<T>::value;


// has width type primary template
template <class T, class U = unsigned int>
struct has_width : std::false_type {};
// Specialization for U = unsigned int
template <class T>
struct has_width <T, decltype((void) T::width, 0u)> : std::true_type {};

// has alignment type primary template
template <class T, class U = unsigned int>
struct has_alignment : std::false_type { };
// Specialization for U = unsigned int
template <class T>
struct has_alignment <T, decltype((void) T::alignment, 0u)> : std::true_type {};

// Patch for missing std::void_t in pre c++17 compilers
template<class... Ts> struct make_void { typedef void type;};
template<class... Ts> using void_t = typename make_void<Ts...>::type;

// Generic template handles types that have no nested ::scalarType member:
template <class, class = void_t<>>
struct has_scalarType : std::false_type { };
// Specialization recognizes types that do have a nested ::scalarType member:
template <class T>
struct has_scalarType<T, void_t<typename T::scalarType>> : std::true_type {};

// If it quacks...
template <class T>
struct is_vector : std::integral_constant
    <bool, has_alignment<T>::value && has_width<T>::value &&
    has_scalarType<T>::value> {};

// Helper c++17 style
// template <class T>
// inline constexpr bool is_vector_v = is_vector<T>::value;


// Generic template handles cases that are not vector type
template <class T, class = void>
struct is_vector_floating_point : std::false_type {};

// Specialized template handles cases that are vector types
template <class T>
struct is_vector_floating_point<T,
    typename std::enable_if<
    is_vector<T>::value>::type
> : std::integral_constant
    <bool, std::is_floating_point<typename T::scalarType>::value> {};

// Helper c++17 style
// template <class T>
// inline constexpr bool is_vector_floating_point_v = is_vector_floating_point<T>::value;

// Generic template handles cases that are not vector type
template <class T, class = void>
struct is_vector_integral : std::false_type {};

// Specialized template handles cases that are vector types
template <class T>
struct is_vector_integral<T,
    typename std::enable_if<
    is_vector<T>::value>::type
> : std::integral_constant
    <bool, std::is_integral<typename T::scalarType>::value> {};

// Helper c++17 style
// template <class T>
// inline constexpr bool is_vector_floating_point_v = is_vector_floating_point<T>::value;


} // namespace tinysimd
#endif
