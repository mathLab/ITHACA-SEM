#pragma once

#include <immintrin.h>
#include <vector>
#include "allocator.hpp"
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


#if defined(__AVX2__) && defined(NEKTAR_ENABLE_SIMD_AVX2)

// forward declaration of concrete types
template<typename T> struct avx2Int8;
template<typename T> struct avx2Long4;
struct avx2Double4;

namespace abi
{

// mapping between abstract types and concrete types
template <> struct avx2<std::int32_t> { using type = avx2Int8<std::int32_t>; };
template <> struct avx2<std::uint32_t> { using type = avx2Int8<std::uint32_t>; };
template <> struct avx2<std::int64_t> { using type = avx2Long4<std::int64_t>; };
template <> struct avx2<std::uint64_t> { using type = avx2Long4<std::uint64_t>; };
template <> struct avx2<double> { using type = avx2Double4; };

} // namespace abi

// concrete types, could add enable if to allow only unsigned long and long...
template<typename T>
struct avx2Int8
{
    static_assert(std::is_integral<T>::value && sizeof(T) == 4,
        "4 bytes Integral required.");

    static constexpr unsigned int width = 8;
    static constexpr unsigned int alignment = 32;

    using scalarType = T;
    using vectorType = __m256i;
    using scalarArray = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline avx2Int8() = default;
    inline avx2Int8(const avx2Int8& rhs) = default;
    inline avx2Int8(const vectorType& rhs) : _data(rhs){}
    inline avx2Int8(const scalarType rhs)
    {
        _data = _mm256_set1_epi32(rhs);
    }
    explicit inline avx2Int8(scalarArray& rhs)
    {
        _data = _mm256_load_si256(reinterpret_cast<vectorType*>(rhs));
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
        _data = _mm256_set1_epi32(rhs);
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

template<typename T>
inline avx2Int8<T> operator+(avx2Int8<T> lhs, avx2Int8<T> rhs)
{
    return _mm256_add_epi32(lhs._data, rhs._data);
}

template<typename T, typename U, typename = typename std::enable_if<
    std::is_arithmetic<U>::value>::type>
inline avx2Int8<T> operator+(avx2Int8<T> lhs, U rhs)
{
    return _mm256_add_epi32(lhs._data, _mm256_set1_epi32(rhs));
}

////////////////////////////////////////////////////////////////////////////////

template<typename T>
struct avx2Long4
{
    static_assert(std::is_integral<T>::value && sizeof(T) == 8,
        "8 bytes Integral required.");

    static constexpr unsigned int width = 4;
    static constexpr unsigned int alignment = 32;

    using scalarType = T;
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
    explicit inline avx2Long4(scalarArray& rhs)
    {
        _data = _mm256_load_si256(reinterpret_cast<vectorType*>(rhs));
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

template<typename T>
inline avx2Long4<T> operator+(avx2Long4<T> lhs, avx2Long4<T> rhs)
{
    return _mm256_add_epi64(lhs._data, rhs._data);
}

template<typename T, typename U, typename = typename std::enable_if<
    std::is_arithmetic<U>::value>::type>
inline avx2Long4<T> operator+(avx2Long4<T> lhs, U rhs)
{
    return _mm256_add_epi64(lhs._data, _mm256_set1_epi64x(rhs));
}

////////////////////////////////////////////////////////////////////////////////

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
    template <typename T>
    inline void gather(scalarType const* p, const sse2Int4<T>& indices)
    {
        _data = _mm256_i32gather_pd(p, indices._data, 8);
    }

    template <typename T>
    inline void gather(scalarType const* p, const avx2Long4<T>& indices)
    {
        _data = _mm256_i64gather_pd(p, indices._data, 8);
    }

    template <typename T>
    inline void scatter(scalarType* out, const sse2Int4<T>& indices) const
    {
        // no scatter intrinsics for AVX2
        alignas(alignment) scalarArray tmp;
        _mm256_store_pd(tmp, _data);

        out[_mm_extract_epi32(indices._data, 0)] = tmp[0];
        out[_mm_extract_epi32(indices._data, 1)] = tmp[1];
        out[_mm_extract_epi32(indices._data, 2)] = tmp[2];
        out[_mm_extract_epi32(indices._data, 3)] = tmp[3];
    }

    template <typename T>
    inline void scatter(scalarType* out, const avx2Long4<T>& indices) const
    {
        // no scatter intrinsics for AVX2
        alignas(alignment) scalarArray tmp;
        _mm256_store_pd(tmp, _data);

        out[_mm256_extract_epi64(indices._data, 0)] = tmp[0];
        out[_mm256_extract_epi64(indices._data, 1)] = tmp[1];
        out[_mm256_extract_epi64(indices._data, 2)] = tmp[2];
        out[_mm256_extract_epi64(indices._data, 3)] = tmp[3];
    }

    // fma
    // this = this + a * b
    inline void fma(const avx2Double4& a, const avx2Double4& b)
    {
        _data = _mm256_fmadd_pd(a._data, b._data, _data);
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

    // unary ops
    inline void operator+=(avx2Double4 rhs)
    {
        _data = _mm256_add_pd(_data, rhs._data);
    }

    inline void operator-=(avx2Double4 rhs)
    {
        _data = _mm256_sub_pd(_data, rhs._data);
    }

    inline void operator*=(avx2Double4 rhs)
    {
        _data = _mm256_mul_pd(_data, rhs._data);
    }

    inline void operator/=(avx2Double4 rhs)
    {
        _data = _mm256_div_pd(_data, rhs._data);
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

inline void load_interleave(
    const double* in,
    size_t dataLen,
    std::vector<avx2Double4, allocator<avx2Double4>> &out)
{
    size_t nBlocks = dataLen / 4;

    alignas(32) size_t tmp[4] = {0, dataLen, 2*dataLen, 3*dataLen};
    using index_t = avx2Long4<size_t>;
    index_t index0(tmp);
    index_t index1 = index0 + 1;
    index_t index2 = index0 + 2;
    index_t index3 = index0 + 3;

    // 4x unrolled loop
    for (size_t i = 0; i < nBlocks; ++i)
    {
        out[4*i + 0].gather(in, index0);
        out[4*i + 1].gather(in, index1);
        out[4*i + 2].gather(in, index2);
        out[4*i + 3].gather(in, index3);
        index0 = index0 + 4;
        index1 = index1 + 4;
        index2 = index2 + 4;
        index3 = index3 + 4;
    }

    // spillover loop
    for (size_t i = 4 * nBlocks; i < dataLen; ++i)
    {
        out[i].gather(in, index0);
        index0 = index0 + 1;
    }
}


inline void deinterleave_store(
    const std::vector<avx2Double4, allocator<avx2Double4>> &in,
    size_t dataLen,
    double *out)
{
#if 0
    double *out0 = out;
    double *out1 = out + dataLen;
    double *out2 = out + 2 * dataLen;
    double *out3 = out + 3 * dataLen;


    for (size_t i = 0; i < dataLen; ++i)
    {
        out0[i] = in[i][0];
        out1[i] = in[i][1];
        out2[i] = in[i][2];
        out3[i] = in[i][3];
    }
#else

    // size_t nBlocks = dataLen / 4;

    alignas(32) size_t tmp[4] = {0, dataLen, 2*dataLen, 3*dataLen};
    using index_t = avx2Long4<size_t>;
    index_t index0(tmp);
    // index_t index1 = index0 + 1;
    // index_t index2 = index0 + 2;
    // index_t index3 = index0 + 3;

    // // 4x unrolled loop
    // for (size_t i = 0; i < nBlocks; ++i)
    // {
    //     in[i].scatter(out, index0);
    //     in[i+1].scatter(out, index1);
    //     in[i+2].scatter(out, index2);
    //     in[i+3].scatter(out, index3);
    //     index0 = index0 + 4;
    //     index1 = index1 + 4;
    //     index2 = index2 + 4;
    //     index3 = index3 + 4;
    // }

    // // spillover loop
    // for (size_t i = 4 * nBlocks; i < dataLen; ++i)
    // {
    //     in[i].scatter(out, index0);
    //     index0 = index0 + 1;
    // }

    for (size_t i = 0; i < dataLen; ++i)
    {
        in[i].scatter(out, index0);
        index0 = index0 + 1;
    }
#endif


}

#endif // defined(__AVX2__)

} // namespace tinysimd
