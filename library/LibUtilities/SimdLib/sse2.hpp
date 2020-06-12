#pragma once

#include <immintrin.h>
#include <cstdint>
#include "traits.hpp"

namespace tinysimd
{

namespace abi
{

template <typename scalarType>
struct sse2
{
    using type = void;
};

} // namespace abi


#if defined(__SSE2__) && defined(NEKTAR_ENABLE_SSE2)

// forward declaration of concrete types
template <typename T>
struct sse2Int4;

namespace abi
{

// mapping between abstract types and concrete types
template <> struct sse2<std::int32_t> { using type = sse2Int4<std::int32_t>; };
template <> struct sse2<std::uint32_t> { using type = sse2Int4<std::uint32_t>; };

} // namespace abi

// concrete types
template <typename T>
struct sse2Int4
{
    static_assert(std::is_integral<T>::value && sizeof(T) == 4,
        "4 bytes Integral required.");

    static constexpr unsigned width = 4;
    static constexpr unsigned alignment = 16;

    using scalarType = T;
    using vectorType = __m128i;
    using scalarArray = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline sse2Int4() = default;
    inline sse2Int4(const sse2Int4& rhs) = default;
    inline sse2Int4(const vectorType& rhs) : _data(rhs){}
    inline sse2Int4(const scalarType rhs)
    {
        _data = _mm_set1_epi32(rhs);
    }

    // store
    inline void store(scalarType* p) const
    {
        _mm_store_si128(reinterpret_cast<vectorType*>(p), _data);
    }

    template
    <
        class flag,
        typename std::enable_if<
            is_requiring_alignment<flag>::value  &&
            !is_streaming<flag>::value, bool
            >::type = 0
    >
    inline void store(scalarType* p, flag) const
    {
        _mm_store_si128(reinterpret_cast<vectorType*>(p), _data);
    }

    template
    <
        class flag,
        typename std::enable_if<
            !is_requiring_alignment<flag>::value, bool
            >::type = 0
    >
    inline void store(scalarType* p, flag) const
    {
        _mm_storeu_si128(reinterpret_cast<vectorType*>(p), _data);
    }

    inline void load(const scalarType* p)
    {
        _data = _mm_load_si128(reinterpret_cast<const vectorType*>(p));
    }

    template
    <
        class flag,
        typename std::enable_if<
            is_requiring_alignment<flag>::value  &&
            !is_streaming<flag>::value, bool
            >::type = 0
    >
    inline void load(const scalarType* p, flag)
    {
        _data = _mm_load_si128(reinterpret_cast<const vectorType*>(p));
    }

    template
    <
        class flag,
        typename std::enable_if<
            !is_requiring_alignment<flag>::value, bool
            >::type = 0
    >
    inline void load(const scalarType* p, flag)
    {
        _data = _mm_loadu_si128(reinterpret_cast<const vectorType*>(p));
    }


    inline void broadcast(const scalarType rhs)
    {
        _data = _mm_set1_epi32(rhs);
    }

    // subscript
    // subscript operators are convienient but expensive
    // should not be used in optimized kernels
    inline scalarType operator[](size_t i) const
    {
        alignas(alignment) scalarArray tmp;
        store(tmp, is_aligned);
        return tmp[i];
    }

    inline scalarType& operator[](size_t i)
    {
        scalarType* tmp = reinterpret_cast<scalarType*>(&_data);
        return tmp[i];
    }


};

#endif // defined(__AVX2__)

} // namespace tinysimd
