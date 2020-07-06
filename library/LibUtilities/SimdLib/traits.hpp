#pragma once

#include <type_traits>

namespace tinysimd
{

// Load tags
static constexpr struct is_aligned_t {} is_aligned;
static constexpr struct is_not_aligned_t {} is_not_aligned;
static constexpr struct is_not_reused_t {} is_not_reused; // streaming, skip cache

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


} // namespace tinysimd