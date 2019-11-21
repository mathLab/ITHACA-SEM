#ifndef VECDATA_HPP
#define VECDATA_HPP

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <boost/align/aligned_allocator.hpp>
#include <immintrin.h>

using namespace Nektar;

template<typename DataType, int vecWidth>
struct VecData
{
};

/////
///// Normal data types
/////
template<typename DataType>
struct VecData<DataType, 1>
{
    using T = VecData<double, 1>;
    static const unsigned int vecWidth = 1;
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
    static const unsigned int vecWidth = 4;
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
using AlignedVector = std::vector<T, boost::alignment::aligned_allocator<T, 64>>;

template<class T, int VW>
AlignedVector<VecData<T, VW>> ToAlignedVector(Array<OneD, T> &input)
{
    size_t nElmt    = input.num_elements();
    int    pad      = nElmt % VW;
    size_t nVecElmt = nElmt + (pad == 0 ? 0 : 1);

    AlignedVector<VecData<T, VW>> ret(nVecElmt);

    T *tmp = &input[0];

    for (int i = 0; i < nElmt / VW; ++i, tmp += VW)
    {
        ret[i] = tmp;
    }

    // Pad everything else out
    if (pad > 0)
    {
        T tmp2[VW];
        for (int i = 0; i < pad; ++i)
        {
            tmp2[i] = input[nElmt / VW + i];
        }
        for (int i = pad; i < VW; ++i)
        {
            tmp2[i] = 0;
        }
        ret[nElmt] = tmp2;
    }

    return ret;
}

#endif
