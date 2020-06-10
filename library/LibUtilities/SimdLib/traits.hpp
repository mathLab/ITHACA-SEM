#pragma once

#include <type_traits>

namespace tinysimd
{

// load tags
static constexpr struct is_aligned_t {} is_aligned;
static constexpr struct is_not_aligned_t {} is_not_aligned;
static constexpr struct is_not_reused_t {} is_not_reused; // streaming, skip cache

template< class T >
struct is_streaming
     : std::integral_constant<
         bool,
         std::is_same<is_not_reused_t, typename std::remove_cv<T>::type>::value
     > {};

template< class T >
struct is_requiring_alignment
     : std::integral_constant<
         bool,
         std::is_same<is_aligned_t, typename std::remove_cv<T>::type>::value ||
         is_streaming<T>::value
     > {};


} // namespace tinysimd

