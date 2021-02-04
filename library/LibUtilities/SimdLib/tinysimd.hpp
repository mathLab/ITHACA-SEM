///////////////////////////////////////////////////////////////////////////////
//
// File: tinysimd.cpp
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
// Description: Light wrapper for automatic selection of available SIMD type.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITES_SIMDLIB_TINYSIMD_H
#define NEKTAR_LIB_LIBUTILITES_SIMDLIB_TINYSIMD_H

#include "scalar.hpp"
#include "avx2.hpp"
#include "avx512.hpp"

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
    using type = typename first_not_void_of
    <
        typename avx512<T>::type,
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
#endif
