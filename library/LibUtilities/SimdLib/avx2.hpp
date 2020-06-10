#pragma once

#include <immintrin.h>
#include "traits.hpp"
#include "sse2.hpp"

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
// struct avx2Int8;
struct avx2Long4;
struct avx2Double4;

namespace abi
{

// mapping between abstract types and concrete types
// template <> struct avx2<int> { using type = avx2Int8; };
template <> struct avx2<long> { using type = avx2Long4; };
template <> struct avx2<size_t> { using type = avx2Long4; };
template <> struct avx2<double> { using type = avx2Double4; };

} // namespace abi

// concrete types
struct avx2Long4
{
    static constexpr unsigned width = 8;
    static constexpr unsigned alignment = 32;

    using scalarType = std::int64_t;
    using vectorType = __m256i;
    using scalarArray = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline avx2Long4() = default;
    inline avx2Long4(const avx2Long4& rhs) = default;
    inline avx2Long4(const vectorType& rhs) : _data(rhs){}
    inline avx2Long4(const scalarType rhs)
    {
        _data = _mm256_set1_epi64x(rhs);
    }

    // store
    inline void store(scalarType* p) const
    {
        _mm256_store_si256(reinterpret_cast<vectorType*>(p), _data);
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
        _mm256_store_si256(reinterpret_cast<vectorType*>(p), _data);
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
        _mm256_storeu_si256(reinterpret_cast<vectorType*>(p), _data);
    }

    inline void load(const scalarType* p)
    {
        _data = _mm256_load_si256(reinterpret_cast<const vectorType*>(p));
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
        _data = _mm256_load_si256(reinterpret_cast<const vectorType*>(p));
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
        _data = _mm256_loadu_si256(reinterpret_cast<const vectorType*>(p));
    }

    inline void broadcast(const scalarType rhs)
    {
        _data = _mm256_set1_epi64x(rhs);
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
    inline void load(scalarType const* a, scalarType const* b,
        scalarType const* c, scalarType const* d)
    {
        __m128d t1, t2, t3, t4;
        __m256d t5;
        t1 = _mm_load_sd(a);                      // SSE2
        t2 = _mm_loadh_pd(t1, b);                 // SSE2
        t3 = _mm_load_sd(c);                      // SSE2
        t4 = _mm_loadh_pd(t3, d);                 // SSE2
        t5 = _mm256_castpd128_pd256(t2);          // cast __m128d -> __m256d
        _data = _mm256_insertf128_pd(t5, t4, 1);
    }

    // broadcast
    inline void broadcast(const scalarType rhs)
    {
        _data = _mm256_set1_pd(rhs);
    }

    //gather
    inline void gather(scalarType const* p, sse2Int4 &indices)
    {
        _data = _mm256_i32gather_pd(p, indices._data, 8);
    }

    inline void gather(scalarType const* p, avx2Long4 &indices)
    {
        _data = _mm256_i64gather_pd(p, indices._data, 8);
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





#endif // defined(__AVX2__)

} // namespace tinysimd
