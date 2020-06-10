#pragma once

#include <immintrin.h>
#include "traits.hpp"

namespace tinysimd
{

namespace abi
{

template <typename scalarType>
struct avx2
{
    using type = void;
};

} // namespace abi


#if defined(__AVX2__)

// forward declaration of concrete types
struct avx2Double4;
// struct avx2Long4;

namespace abi
{

// mapping between abstract types and concrete types
// template <> struct avx2<long> { using type = avx2Long4; };
template <> struct avx2<double> { using type = avx2Double4; };

} // namespace abi

// concrete types
struct avx2Double4
{
    static constexpr unsigned width = 4;
    static constexpr unsigned alignment = 32;

    using scalarType = double;
    using vectorType = __m256d;
    using scalarArray = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline avx2Double4() = default;
    inline avx2Double4(const avx2Double4& rhs) = default;
    inline avx2Double4(const vectorType& rhs) : _data(rhs){}
    inline avx2Double4(const scalarType rhs)
    {
        _data = _mm256_set1_pd(rhs);
    }

    // store
    inline void store(scalarType* p) const
    {
        _mm256_store_pd(p, _data);
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
        _mm256_store_pd(p, _data);
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
        _mm256_storeu_pd(p, _data);
    }

    template
    <
        class flag,
        typename std::enable_if<
            is_streaming<flag>::value, bool
            >::type = 0
    >
    inline void store(scalarType* p, flag) const
    {
        _mm256_stream_pd(p, _data);
    }

    // load packed
    inline void load(const scalarType* p)
    {
        _data = _mm256_load_pd(p);
    }

    template
    <
        class flag,
        typename std::enable_if<
            is_requiring_alignment<flag>::value, bool
            >::type = 0
    >
    inline void load(const scalarType* p, flag)
    {
        _data = _mm256_load_pd(p);
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
        _data = _mm256_loadu_pd(p);
    }

    // load random
    inline void load(const scalarType* a, const scalarType* b,
        const scalarType* c, const scalarType* d)
    {
        __m128d t1, t2, t3, t4;
        __m256d t5;
        t1 = _mm_load_sd(a);                      // SSE
        t2 = _mm_loadh_pd(t1, b);                 // SSE
        t3 = _mm_load_sd(c);                      // SSE
        t4 = _mm_loadh_pd(t3, d);                 // SSE
        t5 = _mm256_castpd128_pd256(t2);          // cast __m128d -> __m256d
        _data = _mm256_insertf128_pd(t5, t4, 1);
    }

    inline void broadcast(const scalarType rhs)
    {
        _data = _mm256_set1_pd(rhs);
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

inline avx2Double4 operator+(avx2Double4 lhs, avx2Double4 rhs)
{
    return _mm256_add_pd(lhs._data, rhs._data);
}

inline avx2Double4 operator-(avx2Double4 lhs, avx2Double4 rhs)
{
    return _mm256_sub_pd(lhs._data, rhs._data);
}

inline avx2Double4 operator*(avx2Double4 lhs, avx2Double4 rhs)
{
    return _mm256_mul_pd(lhs._data, rhs._data);
}

inline avx2Double4 operator/(avx2Double4 lhs, avx2Double4 rhs)
{
    return _mm256_div_pd(lhs._data, rhs._data);
}

inline avx2Double4 sqrt(avx2Double4 in)
{
    return _mm256_sqrt_pd(in._data);
}

inline avx2Double4 abs(avx2Double4 in)
{
    // there is no avx2 _mm256_abs_pd intrinsic
    static const __m256d sign_mask = _mm256_set1_pd(-0.); // -0. = 1 << 63
    return _mm256_andnot_pd(sign_mask, in._data);        // !sign_mask & x
}


// struct avx2_int4: fallback<avx2_int4> {
//     using i32 = std::int32_t;

//     static void copyTo(__m128i v, i32* p) {
//         _mm_storeu_si128((__m128i*)p, v);
//     }
//     static __m128i copyFrom(const i32* p) {
//         return _mm_loadu_si128((const __m128i*)p);
//     }
//     static __m128i broadcast(i32 v) {
//         return _mm_set1_epi32(v);
//     }
//     static __m128i add(__m128i a, __m128i b) {
//         return _mm_add_epi32(a, b);
//     }
//     static __m128i mul(__m128i a, __m128i b) {
//         return _mm_mullo_epi32(a, b);
//     }
//     static __m128i fma(__m128i u, __m128i v, __m128i w) {
//         return add(mul(u, v), w);
//     }
//     static i32 reduce_add(__m128i a) {
//         // Add [a3|a2|a1|a0] to [a2|a3|a0|a1]
//         __m128i b = add(a, _mm_shuffle_epi32(a, 0xb1));
//         // Add [b3|b2|b1|b0] to [b1|b0|b3|b2]
//         __m128i c = add(b, _mm_shuffle_epi32(b, 0x4e));
//         return element(c, 0);
//     }

//     using fallback<avx2_int4>::gather;
//     static __m128i gather(tag<avx2_int4>, const i32* p, __m128i index) {
//         return _mm_i32gather_epi32(p, index, 4);
//     }
// };

#endif // defined(__AVX2__) && defined(__FMA__)

} // namespace tinysimd
