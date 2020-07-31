///////////////////////////////////////////////////////////////////////////////
//
// File: VmathSIMD.hpp
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
// Description: Collection of templated functions for vector mathematics using
// SIMD types, some are functions have optimization tricks such as manual loop
// unrolling.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/SimdLib/tinysimd.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

namespace Vmath
{
namespace SIMD
{
    /// \brief Multiply vector z = x + y
    template<class T,typename = typename std::enable_if
        <
            std::is_floating_point<T>::value
        >::type
    >
    void Vadd(const size_t n, const T *x,  const T *y, T *z)
    {
        using namespace tinysimd;
        using vec_t = simd<T>;

        size_t cnt = n;
        // Vectorized loop unroll 4x
        while (cnt >= 4*vec_t::width)
        {
            // load
            vec_t yChunk0, yChunk1, yChunk2, yChunk3;
            yChunk0.load(y, is_not_aligned);
            yChunk1.load(y + vec_t::width, is_not_aligned);
            yChunk2.load(y + 2*vec_t::width, is_not_aligned);
            yChunk3.load(y + 3*vec_t::width, is_not_aligned);

            vec_t xChunk0, xChunk1, xChunk2, xChunk3;
            xChunk0.load(x, is_not_aligned);
            xChunk1.load(x + vec_t::width, is_not_aligned);
            xChunk2.load(x + 2*vec_t::width, is_not_aligned);
            xChunk3.load(x + 3*vec_t::width, is_not_aligned);

            // z = x + y
            vec_t zChunk0 = xChunk0 + yChunk0;
            vec_t zChunk1 = xChunk1 + yChunk1;
            vec_t zChunk2 = xChunk2 + yChunk2;
            vec_t zChunk3 = xChunk3 + yChunk3;

            // store
            zChunk0.store(z, is_not_aligned);
            zChunk1.store(z + vec_t::width, is_not_aligned);
            zChunk2.store(z + 2*vec_t::width, is_not_aligned);
            zChunk3.store(z + 3*vec_t::width, is_not_aligned);

            // update pointers
            x += 4*vec_t::width;
            y += 4*vec_t::width;
            z += 4*vec_t::width;
            cnt-= 4*vec_t::width;
        }

        // Vectorized loop unroll 2x
        while (cnt >= 2*vec_t::width)
        {
            // load
            vec_t yChunk0, yChunk1;
            yChunk0.load(y, is_not_aligned);
            yChunk1.load(y + vec_t::width, is_not_aligned);

            vec_t xChunk0, xChunk1;
            xChunk0.load(x, is_not_aligned);
            xChunk1.load(x + vec_t::width, is_not_aligned);

            // z = x + y
            vec_t zChunk0 = xChunk0 + yChunk0;
            vec_t zChunk1 = xChunk1 + yChunk1;

            // store
            zChunk0.store(z, is_not_aligned);
            zChunk1.store(z + vec_t::width, is_not_aligned);

            // update pointers
            x += 2*vec_t::width;
            y += 2*vec_t::width;
            z += 2*vec_t::width;
            cnt-= 2*vec_t::width;
        }

        // Vectorized loop
        while (cnt >= vec_t::width)
        {
            // load
            vec_t yChunk;
            yChunk.load(y, is_not_aligned);
            vec_t xChunk;
            xChunk.load(x, is_not_aligned);

            // z = x + y
            vec_t zChunk = xChunk + yChunk;

            // store
            zChunk.store(z, is_not_aligned);

            // update pointers
            x += vec_t::width;
            y += vec_t::width;
            z += vec_t::width;
            cnt-= vec_t::width;
        }

        // spillover loop
        while(cnt)
        {
            // z = x + y;
            *z = (*x) + (*y);
            // update pointers
            ++x;
            ++y;
            ++z;
            --cnt;
        }
    }

    /// \brief Multiply vector z = x * y
    template<class T,typename = typename std::enable_if
        <
            std::is_floating_point<T>::value
        >::type
    >
    void Vmul(const size_t n, const T *x,  const T *y, T *z)
    {
        using namespace tinysimd;
        using vec_t = simd<T>;

        size_t cnt = n;
                // Vectorized loop unroll 4x
        while (cnt >= 4*vec_t::width)
        {
            // load
            vec_t yChunk0, yChunk1, yChunk2, yChunk3;
            yChunk0.load(y, is_not_aligned);
            yChunk1.load(y + vec_t::width, is_not_aligned);
            yChunk2.load(y + 2*vec_t::width, is_not_aligned);
            yChunk3.load(y + 3*vec_t::width, is_not_aligned);

            vec_t xChunk0, xChunk1, xChunk2, xChunk3;
            xChunk0.load(x, is_not_aligned);
            xChunk1.load(x + vec_t::width, is_not_aligned);
            xChunk2.load(x + 2*vec_t::width, is_not_aligned);
            xChunk3.load(x + 3*vec_t::width, is_not_aligned);

            // z = x * y
            vec_t zChunk0 = xChunk0 * yChunk0;
            vec_t zChunk1 = xChunk1 * yChunk1;
            vec_t zChunk2 = xChunk2 * yChunk2;
            vec_t zChunk3 = xChunk3 * yChunk3;

            // store
            zChunk0.store(z, is_not_aligned);
            zChunk1.store(z + vec_t::width, is_not_aligned);
            zChunk2.store(z + 2*vec_t::width, is_not_aligned);
            zChunk3.store(z + 3*vec_t::width, is_not_aligned);

            // update pointers
            x += 4*vec_t::width;
            y += 4*vec_t::width;
            z += 4*vec_t::width;
            cnt-= 4*vec_t::width;
        }

        // Vectorized loop unroll 2x
        while (cnt >= 2*vec_t::width)
        {
            // load
            vec_t yChunk0, yChunk1;
            yChunk0.load(y, is_not_aligned);
            yChunk1.load(y + vec_t::width, is_not_aligned);

            vec_t xChunk0, xChunk1;
            xChunk0.load(x, is_not_aligned);
            xChunk1.load(x + vec_t::width, is_not_aligned);

            // z = x * y
            vec_t zChunk0 = xChunk0 * yChunk0;
            vec_t zChunk1 = xChunk1 * yChunk1;

            // store
            zChunk0.store(z, is_not_aligned);
            zChunk1.store(z + vec_t::width, is_not_aligned);

            // update pointers
            x += 2*vec_t::width;
            y += 2*vec_t::width;
            z += 2*vec_t::width;
            cnt-= 2*vec_t::width;
        }

        // Vectorized loop
        while (cnt >= vec_t::width)
        {
            // load
            vec_t yChunk;
            yChunk.load(y, is_not_aligned);
            vec_t xChunk;
            xChunk.load(x, is_not_aligned);

            // z = x * y
            vec_t zChunk = xChunk * yChunk;

            // store
            zChunk.store(z, is_not_aligned);

            // update pointers
            x += vec_t::width;
            y += vec_t::width;
            z += vec_t::width;
            cnt-= vec_t::width;
        }

        // spillover loop
        while(cnt)
        {
            // z = x * y;
            *z = (*x) * (*y);
            // update pointers
            ++x;
            ++y;
            ++z;
            --cnt;
        }
    }


    /// \brief  vvtvp (vector times vector plus vector): z = w*x + y
    template<class T,typename = typename std::enable_if
        <
            std::is_floating_point<T>::value
        >::type
    >
    void Vvtvp(const size_t n, const T *w,  const T *x,  const T *y, T *z)
    {
        using namespace tinysimd;
        using vec_t = simd<T>;

        size_t cnt = n;
        // Vectorized loop
        while (cnt >= vec_t::width)
        {
            // load
            vec_t wChunk;
            wChunk.load(w, is_not_aligned);
            vec_t yChunk;
            yChunk.load(y, is_not_aligned);
            vec_t xChunk;
            xChunk.load(x, is_not_aligned);

            // z = w * x + y
            vec_t zChunk = wChunk * xChunk + yChunk;

            // store
            zChunk.store(z, is_not_aligned);

            // update pointers
            w += vec_t::width;
            x += vec_t::width;
            y += vec_t::width;
            z += vec_t::width;
            cnt-= vec_t::width;
        }

        // spillover loop
        while(cnt)
        {
            // z = w * x + y;
            *z = (*w) * (*x) + (*y);
            // update pointers
            ++w;
            ++x;
            ++y;
            ++z;
            --cnt;
        }
    }

    /// \brief  vvtvm (vector times vector plus vector): z = w*x - y
    template<class T,typename = typename std::enable_if
        <
            std::is_floating_point<T>::value
        >::type
    >
    void Vvtvm(const size_t n, const T *w,  const T *x,  const T *y, T *z)
    {
        using namespace tinysimd;
        using vec_t = simd<T>;

        size_t cnt = n;
        // Vectorized loop
        while (cnt >= vec_t::width)
        {
            // load
            vec_t wChunk;
            wChunk.load(w, is_not_aligned);
            vec_t yChunk;
            yChunk.load(y, is_not_aligned);
            vec_t xChunk;
            xChunk.load(x, is_not_aligned);

            // z = w * x - y
            vec_t zChunk = wChunk * xChunk - yChunk;

            // store
            zChunk.store(z, is_not_aligned);

            // update pointers
            w += vec_t::width;
            x += vec_t::width;
            y += vec_t::width;
            z += vec_t::width;
            cnt-= vec_t::width;
        }

        // spillover loop
        while(cnt)
        {
            // z = w * x - y;
            *z = (*w) * (*x) - (*y);
            // update pointers
            ++w;
            ++x;
            ++y;
            ++z;
            --cnt;
        }
    }

    /// \brief  vvtvvtp (vector times vector plus vector times vector):
    // z = v*w + x*y
    template<class T,typename = typename std::enable_if
        <
            std::is_floating_point<T>::value
        >::type
    >
    inline void Vvtvvtp (const size_t n, const T* v, const T* w, const T* x,
        const T* y, T* z)
    {
        using namespace tinysimd;
        using vec_t = simd<T>;

        size_t cnt = n;
        // Vectorized loop
        while (cnt >= vec_t::width)
        {
            // load
            vec_t vChunk;
            vChunk.load(v, is_not_aligned);
            vec_t wChunk;
            wChunk.load(w, is_not_aligned);
            vec_t yChunk;
            yChunk.load(y, is_not_aligned);
            vec_t xChunk;
            xChunk.load(x, is_not_aligned);

            // z = v * w + x * y;
            vec_t z1Chunk = vChunk * wChunk;
            vec_t z2Chunk = xChunk * yChunk;
            vec_t zChunk = z1Chunk + z2Chunk;

            // store
            zChunk.store(z, is_not_aligned);

            // update pointers
            v += vec_t::width;
            w += vec_t::width;
            x += vec_t::width;
            y += vec_t::width;
            z += vec_t::width;
            cnt-= vec_t::width;
        }

        // spillover loop
        while(cnt)
        {
            // z = v * w + x * y;
            T z1 = (*v) * (*w);
            T z2 = (*x) * (*y);
            *z = z1 + z2;
            // update pointers
            ++v;
            ++w;
            ++x;
            ++y;
            ++z;
            --cnt;
        }
    }

    /// \brief Gather vector z[i] = x[y[i]]
    template<class T, class I, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value &&
            std::is_integral<I>::value
        >::type
    >
    void Gathr(const I n, const T* x,  const I* y, T* z)
    {
        using namespace tinysimd;
        using vec_t = simd<T>;
        using vec_t_i = simd<I>;

        I cnt = n;
        // Unroll 4x Vectorized loop
        while (cnt >= 4*vec_t::width)
        {
            // load index
            vec_t_i yChunk0, yChunk1, yChunk2, yChunk3;
            yChunk0.load(y, is_not_aligned);
            yChunk1.load(y + vec_t_i::width, is_not_aligned);
            yChunk2.load(y + 2*vec_t_i::width, is_not_aligned);
            yChunk3.load(y + 3*vec_t_i::width, is_not_aligned);

            // z = x[y[i]]
            vec_t zChunk0, zChunk1, zChunk2, zChunk3;
            zChunk0.gather(x, yChunk0);
            zChunk1.gather(x, yChunk1);
            zChunk2.gather(x, yChunk2);
            zChunk3.gather(x, yChunk3);

            // store
            zChunk0.store(z, is_not_aligned);
            zChunk1.store(z + vec_t_i::width, is_not_aligned);
            zChunk2.store(z + 2*vec_t_i::width, is_not_aligned);
            zChunk3.store(z + 3*vec_t_i::width, is_not_aligned);

            // update pointers
            y += 4*vec_t_i::width;
            z += 4*vec_t::width;
            cnt-= 4*vec_t::width;
        }

        // Unroll 2x Vectorized loop
        while (cnt >= 2*vec_t::width)
        {
            // load index
            vec_t_i yChunk0, yChunk1;
            yChunk0.load(y, is_not_aligned);
            yChunk1.load(y + vec_t_i::width, is_not_aligned);

            // z = x[y[i]]
            vec_t zChunk0, zChunk1;
            zChunk0.gather(x, yChunk0);
            zChunk1.gather(x, yChunk1);

            // store
            zChunk0.store(z, is_not_aligned);
            zChunk1.store(z + vec_t_i::width, is_not_aligned);

            // update pointers
            y += 2*vec_t_i::width;
            z += 2*vec_t::width;
            cnt-= 2*vec_t::width;
        }

        // Vectorized loop
        while (cnt >= vec_t::width)
        {
            // load index
            vec_t_i yChunk;
            yChunk.load(y, is_not_aligned);

            // z = x[y[i]]
            vec_t zChunk;
            zChunk.gather(x, yChunk);

            // store
            zChunk.store(z, is_not_aligned);

            // update pointers
            y += vec_t_i::width;
            z += vec_t::width;
            cnt-= vec_t::width;
        }

        // spillover loop
        while(cnt)
        {
            // z = x[y[i]]
            *z = *(x + *y);
            // update pointers
            ++y;
            ++z;
            --cnt;
        }
    }





} // namespace SIMD
} // namespace Vmath

