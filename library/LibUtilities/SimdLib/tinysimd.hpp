#pragma once

#include "avx2.hpp"
#include "scalar.hpp"

namespace tinysimd
{

template <typename...>
struct first_not_void_of { using type = void; };

template <typename... Rest>
struct first_not_void_of<void, Rest...>
{
    using type = typename first_not_void_of<Rest...>::type;
};

template <typename T, typename... Rest>
struct first_not_void_of<T, Rest...>
{
    using type = T;
};

namespace abi
{

// pick the most specialiazed available match out of available ABIs.
template <typename T> struct default_abi
{
    using type = typename first_not_void_of<
        typename avx2<T>::type,
        typename sse2<T>::type,
        typename scalar<T>::type
    >::type;

    static_assert(!std::is_void<type>::value, "unsupported SIMD type");
};

} // namespace abi


// light wrapper for default types
template <typename ScalarType,
    template <typename> class abi = abi::default_abi>
using simd = typename abi<ScalarType>::type;



} // namespace tinysimd

