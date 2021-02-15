///////////////////////////////////////////////////////////////////////////////
//
// File: scalar.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Scalar type used when a vector type is needed, but no SIMD
// extension is available.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITES_SIMDLIB_SCALAR_H
#define NEKTAR_LIB_LIBUTILITES_SIMDLIB_SCALAR_H

#include <vector>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <cstdint>
#include "allocator.hpp"
#include "traits.hpp"

namespace tinysimd
{

namespace abi
{

template <typename scalarType>
struct scalar
{
    using type = void;
};

} // namespace abi

// forward declaration of concrete types
// makes default type available for all arithmetic types
template<typename T, typename = typename std::enable_if<
    std::is_arithmetic<T>::value>::type>
struct scalarT;
struct scalarMask;

namespace abi
{

// mapping between abstract types and concrete types
template <> struct scalar<double>{ using type = scalarT<double>; };
template <> struct scalar<std::int64_t>{ using type = scalarT<std::int64_t>; };
template <> struct scalar<std::uint64_t>{ using type = scalarT<std::uint64_t>; };
template <> struct scalar<std::int32_t>{ using type = scalarT<std::int32_t>; };
template <> struct scalar<std::uint32_t>{ using type = scalarT<std::uint32_t>; };
template <> struct scalar<bool>{ using type = scalarMask; };
template <> struct scalar<size_t>{ using type = scalarT<size_t>; };

} // namespace abi

// concrete types
template<typename T, typename>
struct scalarT
{
    static constexpr unsigned int width = 1;
    static constexpr unsigned int alignment = sizeof(T);

    using scalarType = T;
    using vectorType = scalarType;
    using scalarArray = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline scalarT() = default;
    inline scalarT(const scalarT& rhs) = default;
    inline scalarT(const vectorType& rhs) : _data(rhs){}

    // store
    inline void store(scalarType* p) const
    {
        *p = _data;
    }

    template<class flag>
    inline void store(scalarType* p, flag) const
    {
        *p = _data;
    }

    // load
    inline void load(const scalarType* p)
    {
        _data = *p;
    }

    template<class flag>
    inline void load(const scalarType* p, flag)
    {
        _data = *p;
    }

    inline void broadcast(const scalarType rhs)
    {
        _data = rhs;
    }

    template<typename U, typename = typename std::enable_if<
        std::is_integral<U>::value>::type>
    inline void gather(const scalarType* p, const scalarT<U>& indices)
    {
        _data = *(p + indices._data);
    }

    template<typename U, typename = typename std::enable_if<
        std::is_integral<U>::value>::type>
    inline void scatter(scalarType* p, const scalarT<U>& indices) const
    {
        p += indices._data;
        *p = _data;
    }

    // fma
    // this = this + a * b
    inline void fma(const scalarT<T>& a, const scalarT<T>& b)
    {
        _data += a._data * b._data;
    }

    // subscript
    inline scalarType operator[](size_t) const
    {
        return _data;
    }

    inline scalarType& operator[](size_t)
    {
        return _data;
    }

    // unary ops
    inline void operator+=(scalarT<T> rhs)
    {
        _data += rhs._data;
    }

    inline void operator-=(scalarT<T> rhs)
    {
        _data -= rhs._data;
    }

    inline void operator*=(scalarT<T> rhs)
    {
        _data *= rhs._data;
    }

    inline void operator/=(scalarT<T> rhs)
    {
        _data /= rhs._data;
    }

};

template<typename T>
inline scalarT<T> operator+(scalarT<T> lhs, scalarT<T> rhs)
{
    return lhs._data + rhs._data;
}
template<typename T, typename U, typename = typename std::enable_if<
    std::is_arithmetic<U>::value>::type>
inline scalarT<T> operator+(U lhs, scalarT<T> rhs)
{
    return lhs + rhs._data;
}
template<typename T, typename U, typename = typename std::enable_if<
    std::is_arithmetic<U>::value>::type>
inline scalarT<T> operator+(scalarT<T> lhs, U rhs)
{
    return lhs._data + rhs;
}


template<typename T>
inline scalarT<T> operator-(scalarT<T> lhs, scalarT<T> rhs)
{
    return lhs._data - rhs._data;
}
template<typename T, typename U, typename = typename std::enable_if<
    std::is_arithmetic<U>::value>::type>
inline scalarT<T> operator-(U lhs, scalarT<T> rhs)
{
    return lhs - rhs._data;
}
template<typename T, typename U, typename = typename std::enable_if<
    std::is_arithmetic<U>::value>::type>
inline scalarT<T> operator-(scalarT<T> lhs, U rhs)
{
    return lhs._data - rhs;
}


template<typename T>
inline scalarT<T> operator*(scalarT<T> lhs, scalarT<T> rhs)
{
    return lhs._data * rhs._data;
}
template<typename T, typename U, typename = typename std::enable_if<
    std::is_arithmetic<U>::value>::type>
inline scalarT<T> operator*(U lhs, scalarT<T> rhs)
{
    return lhs * rhs._data;
}
template<typename T, typename U, typename = typename std::enable_if<
    std::is_arithmetic<U>::value>::type>
inline scalarT<T> operator*(scalarT<T> lhs, U rhs)
{
    return lhs._data * rhs;
}


template<typename T>
inline scalarT<T> operator/(scalarT<T> lhs, scalarT<T> rhs)
{
    return lhs._data / rhs._data;
}
template<typename T, typename U, typename = typename std::enable_if<
    std::is_arithmetic<U>::value>::type>
inline scalarT<T> operator/(U lhs, scalarT<T> rhs)
{
    return lhs / rhs._data;
}
template<typename T, typename U, typename = typename std::enable_if<
    std::is_arithmetic<U>::value>::type>
inline scalarT<T> operator/(scalarT<T> lhs, U rhs)
{
    return lhs._data / rhs;
}

template<typename T>
inline scalarT<T> sqrt(scalarT<T> in)
{
    return std::sqrt(in._data);
}
template<typename T>
inline scalarT<T> abs(scalarT<T> in)
{
    return std::abs(in._data);
}

template<typename T>
inline scalarT<T> log(scalarT<T> in)
{
    return std::log(in._data);
}

template<typename T>
inline void load_interleave(
    const T* in,
    size_t dataLen,
    std::vector<scalarT<T>, allocator<scalarT<T>>> &out)
{
    for (size_t i = 0; i < dataLen; ++i)
    {
        out[i] = in[i];
    }
}

template<typename T>
inline void deinterleave_store(
    const std::vector<scalarT<T>, allocator<scalarT<T>>> &in,
    size_t dataLen,
    T *out)
{
    for (size_t i = 0; i < dataLen; ++i)
    {
        out[i] = in[i]._data;
    }
}

////////////////////////////////////////////////////////////////////////////////

// mask type
// mask is a int type that uses boolean promotion
//
// VERY LIMITED SUPPORT...just enough to make cubic eos work...
//
struct scalarMask : scalarT<std::uint64_t>
{
    // bring in ctors
    using scalarT::scalarT;

    static constexpr scalarType true_v = true;
    static constexpr scalarType false_v = false;
};

inline scalarMask operator>(scalarT<double> lhs, scalarT<double> rhs)
{

    return lhs._data > rhs._data;
}

inline bool operator&&(scalarMask lhs, bool rhs)
{
    return lhs._data && rhs;
}

} // namespace tinysimd
#endif
