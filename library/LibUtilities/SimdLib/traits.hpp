#pragma once

#include <type_traits>

namespace tinysimd
{

// load tags
static constexpr struct is_aligned_t {} is_aligned;
static constexpr struct is_not_aligned_t {} is_not_aligned;
static constexpr struct is_not_reused_t {} is_not_reused; // streaming, skip cache

// check load tags
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

// helper c++17 style
// template< class T >
// inline constexpr bool is_streaming_v = is_streaming<T>::value;
// template< class T >
// inline constexpr bool is_requiring_alignment_v = is_requiring_alignment<T>::value;


// has width type primary template
template <typename T, typename U = unsigned int>
struct has_width : std::false_type { };
// Specialization for U = int
template <typename T>
struct has_width <T, decltype((void) T::width, 0u)> : std::true_type {};

// has alignment type primary template
template <typename T, typename U = unsigned int>
struct has_alignment : std::false_type { };
// Specialization for U = int
template <typename T>
struct has_alignment <T, decltype((void) T::alignment, 0u)> : std::true_type {};


// if it quacks...
template <class T>
struct is_vector : std::integral_constant
    <bool, has_alignment<T>::value && has_width<T>::value> {};

// helper c++17 style
// template< class T >
// inline constexpr bool is_vector_v = is_vector<T>::value;


} // namespace tinysimd