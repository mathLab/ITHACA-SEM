#ifndef VECDATA_HPP
#define VECDATA_HPP

#include <boost/align/aligned_allocator.hpp>
#include <vector>
#include <iostream>
#include <immintrin.h>

namespace Nektar {
namespace AVX {

#if defined(__AVX2__)
/// Number of bits in a AVX2 vector
constexpr int SIMD_WIDTH_BITS = 256;
#endif
#if defined(__AVX512F__)
/// Number of bits in a AVX512 vector
constexpr int SIMD_WIDTH_BITS = 512;
#endif

/// Number of bytes in a AVX vector
constexpr int SIMD_WIDTH_BYTES = SIMD_WIDTH_BITS / 8;
/// Number of elements in a AVX vector for NekDoubles
constexpr int SIMD_WIDTH_SIZE = SIMD_WIDTH_BYTES / sizeof(double);


template<typename DataType, int vecWidth>
struct VecData
{
};

template<class T>
using AlignedVector =
    std::vector<T, boost::alignment::aligned_allocator<T, SIMD_WIDTH_BYTES>>;

/////
///// Normal data types
/////
template<typename DataType>
struct VecData<DataType, 1>
{
    using T = VecData<double, 1>;
    static const size_t vecWidth = 1;
    DataType m_data;

    inline VecData() = default;
    inline VecData(const DataType &rhs)
    {
        m_data = rhs;
    }
    inline VecData(const DataType *rhs)
    {
        m_data = *rhs;
    }

    inline T &operator=(const DataType &data)
    {
        m_data = data;
        return *this;
    }

    inline T &operator=(const DataType *data)
    {
        m_data = *data;
        return *this;
    }

    inline void fma(const T &a, const T &b)
    {
        m_data += a.m_data * b.m_data;
    }

    inline void store(DataType *out)
    {
        *out = m_data;
    }

    inline void store_nts(double *out)
    {
        *out = m_data;
    }

    void scatter(DataType *out, VecData<int,1> index)
    {
        out[index] = m_data;
    }

    inline static void load_interleave(
        DataType *in,
        size_t dataLen,
        AlignedVector<VecData<DataType, 1>> &out)
    {
        // Nothing to do -- straightforward copy.
        for (size_t i = 0; i < dataLen; ++i)
        {
            out[i] = in[i];
        }
    }

};

template<typename DataType>
inline VecData<DataType, 1> operator*(VecData<DataType, 1> lhs, VecData<DataType, 1> rhs)
{
    return lhs.m_data * rhs.m_data;
}

template<typename DataType>
inline VecData<DataType, 1> operator+(VecData<DataType, 1> lhs, VecData<DataType, 1> rhs)
{
    return lhs.m_data + rhs.m_data;
}

template<typename DataType>
inline VecData<DataType, 1> operator-(VecData<DataType, 1> lhs, VecData<DataType, 1> rhs)
{
    return lhs.m_data - rhs.m_data;
}

template<typename DataType>
inline VecData<DataType, 1> operator/(VecData<DataType, 1> lhs, VecData<DataType, 1> rhs)
{
    return lhs.m_data / rhs.m_data;
}

template<typename DataType>
inline VecData<DataType, 1> gather(const DataType* data, VecData<int, 1> &index)
{
    return data[index];
}

template<typename DataType>
inline std::ostream &operator<<(std::ostream &os, VecData<DataType, 1> const &vec) {
    return os << vec.m_data;
}
#if defined(__AVX2__)
/////
///// Specialisation: int(32-bit) + AVX
/////
template<>
struct VecData<int, 4>
{
    using T = VecData<int,4>;
    static const size_t vecWidth = 4;
    __m128i m_data;

    inline VecData() = default;
    inline VecData(const T &rhs) = default;
    inline VecData(const int &rhs)
    {
        m_data = _mm_set1_epi32(rhs);
    }
    inline VecData(const int *rhs)
    {
        m_data = _mm_load_si128((__m128i*)rhs);
    }
    inline VecData(const __m128i &rhs) : m_data(rhs)
    {
    }
    inline VecData& operator=(const __m128i &rhs)
    {
        m_data = rhs;
        return *this;
    }
    inline T &operator=(const int &data)
    {
        m_data = _mm_set1_epi32(data);
        return *this;
    }
    inline T &operator=(const int *data)
    {
        m_data = _mm_load_si128((__m128i*)data);
        return *this;
    }

};

inline std::ostream &operator<<(std::ostream &os, VecData<int, 4> const &vec) {
    int d[4];
    _mm_store_si128((__m128i*)d, vec.m_data);
    return os << d[0] << ", " << d[1] << ", " << d[2] << ", " << d[3] << ";";
}

/////
///// Specialisation: double + AVX2 (including FMA)
/////
template<>
struct VecData<double, 4>
{
    using T = VecData<double, 4>;
    static const unsigned int vecWidth = 4;
    __m256d m_data;

    inline VecData() = default;
    inline VecData(const T &rhs) = default;
    inline VecData(const double &rhs)
    {
        m_data = _mm256_set1_pd(rhs);
    }
    inline VecData(const double *rhs)
    {
        m_data = _mm256_load_pd(rhs);
    }
    inline VecData(const __m256d &rhs) : m_data(rhs)
    {
    }

    inline VecData& operator=(const __m256d &rhs)
    {
        m_data = rhs;
        return *this;
    }
    inline T &operator=(const double &data)
    {
        m_data = _mm256_set1_pd(data);
        return *this;
    }

    inline T &operator=(const double *data)
    {
        m_data = _mm256_load_pd(data);
        return *this;
    }

    inline void fma(const T &a, const T &b)
    {
        m_data = _mm256_fmadd_pd(a.m_data, b.m_data, m_data);
    }

    inline void store(double *out)
    {
        _mm256_store_pd(out, m_data);
    }

    inline void store_nts(double *out)
    {
        _mm256_stream_pd(out, m_data);
    }

    void scatter(double* out, VecData<int,4> indices)
    {
        //No scatter inrins for AVX2, so we just do it manually

        double d[4]; //Couldn't find an extract for 1 double
        _mm256_store_pd(d, m_data);

        out[_mm_extract_epi32(indices.m_data,0)] = d[0];
        out[_mm_extract_epi32(indices.m_data,1)] = d[1];
        out[_mm_extract_epi32(indices.m_data,2)] = d[2];
        out[_mm_extract_epi32(indices.m_data,3)] = d[3];
    }


    inline static void load_interleave(
        const double *in,
        size_t dataLen,
        AlignedVector<VecData<double, 4>> &out)
    {
        size_t nBlocks = dataLen / 4;

        const double *d0 = in;
        const double *d1 = in + dataLen;
        const double *d2 = in + 2 * dataLen;
        const double *d3 = in + 3 * dataLen;

        /*
        for (size_t i = 0; i < dataLen; ++i)
        {
            out[i].m_data[0] = d0[i];
            out[i].m_data[1] = d1[i];
            out[i].m_data[2] = d2[i];
            out[i].m_data[3] = d3[i];
        }
        return;
        */

        for (size_t i = 0; i < nBlocks; ++i)
        {
#if 0
            // Load 4 doubles from each of d0..d3 to give 2d array:
            //
            // [ d00 d01 d02 d03 ] => ymm0
            // [ d10 d11 d12 d13 ] => ymm1
            // [ d20 d21 d22 d23 ] => ymm2
            // [ d30 d31 d32 d33 ] => ymm3
            //
            // Now want to now transpose this and store in contiguous locations
            // in out[] array.
            __m256d ymm0 = _mm256_load_pd(d0 + 4*i);
            __m256d ymm1 = _mm256_load_pd(d1 + 4*i);
            __m256d ymm2 = _mm256_load_pd(d2 + 4*i);
            __m256d ymm3 = _mm256_load_pd(d3 + 4*i);

            // Perform a shuffle so that e.g.. The below uses vperm2f128 which
            // is not optimal since it is only executed on port 5 of
            // processors. See below for an implementation that uses
            // _mm256_insertf128 instead.
            __m256d ymm4 = _mm256_permute2f128_pd(ymm0, ymm2, 0x20); // ymm4 <= [d00 d01 d20 d21]
            __m256d ymm5 = _mm256_permute2f128_pd(ymm1, ymm3, 0x20); // ymm5 <= [d10 d11 d30 d31]
            __m256d ymm6 = _mm256_permute2f128_pd(ymm0, ymm2, 0x31); // ymm6 <= [d02 d03 d22 d23]
            __m256d ymm7 = _mm256_permute2f128_pd(ymm1, ymm3, 0x31); // ymm7 <= [d12 d13 d32 d33]
#else
            __m256d ymm4 = _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm_loadu_pd(d0)), _mm_loadu_pd(d2), 1);
            __m256d ymm5 = _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm_loadu_pd(d1)), _mm_loadu_pd(d3), 1);
            __m256d ymm6 = _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm_loadu_pd(d0+2)), _mm_loadu_pd(d2+2), 1);
            __m256d ymm7 = _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm_loadu_pd(d1+2)), _mm_loadu_pd(d3+2), 1);
#endif

            // Simultaneous unpack/interleave
            out[4*i + 0] = _mm256_unpacklo_pd(ymm4, ymm5); // out <= [d00 d10 d20 d30]
            out[4*i + 1] = _mm256_unpackhi_pd(ymm4, ymm5); // out <= [d01 d11 d21 d31]
            out[4*i + 2] = _mm256_unpacklo_pd(ymm6, ymm7); // out <= [d02 d12 d22 d32]
            out[4*i + 3] = _mm256_unpackhi_pd(ymm6, ymm7); // out <= [d03 d13 d23 d33]

            d0 += 4;
            d1 += 4;
            d2 += 4;
            d3 += 4;
        }

        // Handle leftover data.
        for (size_t i = 4 * nBlocks; i < dataLen; ++i)
        {
            out[i].m_data[0] = *d0;
            out[i].m_data[1] = *d1;
            out[i].m_data[2] = *d2;
            out[i].m_data[3] = *d3;
            ++d0;
            ++d1;
            ++d2;
            ++d3;
        }
    }

    inline static void deinterleave_store(
        const AlignedVector<VecData<double, 4>> &in,
        size_t dataLen,
        double *out)
    {
        // unsigned int nBlocks = dataLen / 4;

        double *out0 = out;
        double *out1 = out + dataLen;
        double *out2 = out + 2 * dataLen;
        double *out3 = out + 3 * dataLen;

        /*
        for (size_t i = 0; i < nBlocks; ++i)
        {
            __m256d ymm4 = _mm256_permute2f128_pd(in[4*i+0].m_data, in[4*i+2].m_data, 0x20); // ymm4 <= [d00 d01 d20 d21]
            __m256d ymm5 = _mm256_permute2f128_pd(in[4*i+1].m_data, in[4*i+3].m_data, 0x20); // ymm5 <= [d10 d11 d30 d31]
            __m256d ymm6 = _mm256_permute2f128_pd(in[4*i+0].m_data, in[4*i+2].m_data, 0x31); // ymm6 <= [d02 d03 d22 d23]
            __m256d ymm7 = _mm256_permute2f128_pd(in[4*i+1].m_data, in[4*i+3].m_data, 0x31); // ymm7 <= [d12 d13 d32 d33]

            _mm256_store_pd(out0, _mm256_unpacklo_pd(ymm4, ymm5));
            _mm256_store_pd(out1, _mm256_unpackhi_pd(ymm4, ymm5));
            _mm256_store_pd(out2, _mm256_unpacklo_pd(ymm6, ymm7));
            _mm256_store_pd(out3, _mm256_unpackhi_pd(ymm6, ymm7));

            out0 += 4;
            out1 += 4;
            out2 += 4;
            out3 += 4;
        }

        for (size_t i = 4 * nBlocks; i < dataLen; ++i)
        {
            *out0 = in[i].m_data[0];
            *out1 = in[i].m_data[1];
            *out2 = in[i].m_data[2];
            *out3 = in[i].m_data[3];
            ++out0;
            ++out1;
            ++out2;
            ++out3;
        }
        */
        for (size_t i = 0; i < dataLen; ++i)
        {
            out0[i] = in[i].m_data[0];
            out1[i] = in[i].m_data[1];
            out2[i] = in[i].m_data[2];
            out3[i] = in[i].m_data[3];
        }
    }
};

inline VecData<double, 4> operator*(VecData<double, 4> lhs, VecData<double, 4> rhs)
{
    return _mm256_mul_pd(lhs.m_data, rhs.m_data);
}

inline VecData<double, 4> operator+(VecData<double, 4> lhs, VecData<double, 4> rhs)
{
    return _mm256_add_pd(lhs.m_data, rhs.m_data);
}

inline VecData<double, 4> operator-(VecData<double, 4> lhs, VecData<double, 4> rhs)
{
    return _mm256_sub_pd(lhs.m_data, rhs.m_data);
}

inline VecData<double, 4> operator/(VecData<double, 4> lhs, VecData<double, 4> rhs)
{
    return _mm256_div_pd(lhs.m_data, rhs.m_data);
}

inline std::ostream &operator<<(std::ostream &os, VecData<double, 4> const &vec) {
    double d[4];
    _mm256_store_pd(d, vec.m_data);
    return os << d[0] << ", " << d[1] << ", " << d[2] << ", " << d[3] << ";";
}

//Gather/Scatter operations
inline VecData<double, 4> gather(const double* data, VecData<int, 4> &indices)
{
    return VecData<double,4>( _mm256_i32gather_pd(data, indices.m_data, 8) );
}

// sqrt
inline VecData<double, 4> sqrt(VecData<double, 4> in)
{
    return _mm256_sqrt_pd(in.m_data);
}

// abs
inline VecData<double, 4> abs(VecData<double, 4> in)
{
    // there is no avx2 _mm256_abs_pd intrinsic
    static const __m256d sign_mask = _mm256_set1_pd(-0.); // -0. = 1 << 63
    return _mm256_andnot_pd(sign_mask, in.m_data);        // !sign_mask & x
}


#endif

#if defined(__AVX512F__)
/////
///// Specialisation: int(32-bit) + AVX
/////
template<>
struct VecData<int, 8>
{
    using T = VecData<int,8>;
    static const unsigned int vecWidth = 8;
    __m256i m_data;

    inline VecData() = default;
    inline VecData(const T &rhs) = default;
    inline VecData(const int &rhs)
    {
        m_data = _mm256_set1_epi32(rhs);
    }
    inline VecData(const int *rhs)
    {
        m_data = _mm256_load_si256((__m256i*)rhs);
    }
    inline VecData(const __m256i &rhs) : m_data(rhs)
    {
    }
    inline VecData& operator=(const __m256i &rhs)
    {
        m_data = rhs;
        return *this;
    }
    inline T &operator=(const int &data)
    {
        m_data = _mm256_set1_epi32(data);
        return *this;
    }
    inline T &operator=(const int *data)
    {
        m_data = _mm256_load_si256((__m256i*)data);
        return *this;
    }

};


/////
///// Specialisation: double + AVX-512
/////
template<>
struct VecData<double, 8>
{
    using T = VecData<double, 8>;
    static const unsigned int vecWidth = 8;
    __m512d m_data;

    inline VecData() = default;
    inline VecData(const T &rhs) = default;
    inline VecData(const double &rhs)
    {
        m_data = _mm512_set1_pd(rhs);
    }
    inline VecData(const double *rhs)
    {
        m_data = _mm512_load_pd(rhs);
    }
    inline VecData(const __m512d &rhs) : m_data(rhs)
    {
    }

    inline VecData& operator=(const __m512d &rhs)
    {
        m_data = rhs;
        return *this;
    }
    inline T &operator=(const double &data)
    {
        m_data = _mm512_set1_pd(data);
        return *this;
    }

    inline T &operator=(const double *data)
    {
        m_data = _mm512_load_pd(data);
        return *this;
    }

    inline void fma(const T &a, const T &b)
    {
        m_data = _mm512_fmadd_pd(a.m_data, b.m_data, m_data);
    }

    inline void store(double *out)
    {
        _mm512_store_pd(out, m_data);
    }

    inline void store_nts(double *out)
    {
        _mm512_stream_pd(out, m_data);
    }

    void scatter(double* out, VecData<int,8> indices)
    {
        _mm512_i32scatter_pd(out, indices.m_data, 8);
    }

};

inline VecData<double, 8> operator*(VecData<double, 8> lhs, VecData<double, 8> rhs)
{
    return _mm512_mul_pd(lhs.m_data, rhs.m_data);
}

inline VecData<double, 8> operator+(VecData<double, 8> lhs, VecData<double, 8> rhs)
{
    return _mm512_add_pd(lhs.m_data, rhs.m_data);
}

inline VecData<double, 8> operator-(VecData<double, 8> lhs, VecData<double, 8> rhs)
{
    return _mm512_sub_pd(lhs.m_data, rhs.m_data);
}

inline VecData<double, 8> operator/(VecData<double, 8> lhs, VecData<double, 8> rhs)
{
    return _mm512_div_pd(lhs.m_data, rhs.m_data);
}

inline std::ostream &operator<<(std::ostream &os, VecData<double, 8> const &vec) {
    double d[8];
    _mm512_store_pd(d, vec.m_data);
    return os << d[0] << ", " << d[1] << ", " << d[2] << ", " << d[3] << ", "
              << d[4] << ", " << d[5] << ", " << d[6] << ", " << d[7] << ";";
}

inline VecData<double, 8> gather(const double* data, VecData<int, 8> &indices)
{
    return VecData<double,8>( _mm512_i32gather_pd(indices.m_data, data, 8) );
}

#endif

/////
///// Aligned vector
/////

template<class T>
using AlignedVector =
    std::vector<T, boost::alignment::aligned_allocator<T, SIMD_WIDTH_BYTES>>;

// template<class T, int VW>
// AlignedVector<VecData<T, VW>> ToAlignedVector(Array<OneD, T> &input)
// {
//     size_t nElmt    = input.num_elements();
//     int    pad      = nElmt % VW;
//     size_t nVecElmt = nElmt + (pad == 0 ? 0 : 1);

//     AlignedVector<VecData<T, VW>> ret(nVecElmt);

//     T *tmp = &input[0];

//     for (int i = 0; i < nElmt / VW; ++i, tmp += VW)
//     {
//         ret[i] = tmp;
//     }

//     // Pad everything else out
//     if (pad > 0)
//     {
//         T tmp2[VW];
//         for (int i = 0; i < pad; ++i)
//         {
//             tmp2[i] = input[nElmt / VW + i];
//         }
//         for (int i = pad; i < VW; ++i)
//         {
//             tmp2[i] = 0;
//         }
//         ret[nElmt] = tmp2;
//     }

//     return ret;
// }

} // namespace AVX
} // namespace Nektar

#endif // header guard