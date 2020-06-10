///////////////////////////////////////////////////////////////////////////////
//
// File: TestVmathSIMD.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VmathSIMD.hpp>

#include <boost/core/ignore_unused.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
namespace VmathSIMDUnitTests
{
    BOOST_AUTO_TEST_CASE(TestVvtvp)
    {
        using dataType = double;
        constexpr size_t n = 11;
        alignas(tinysimd::simd<dataType>::alignment)
            std::array<dataType, n> w, x, y, z;
        dataType epsilon = std::numeric_limits<dataType>::epsilon();

        // init
        for (size_t i = 0; i < n; ++i)
        {
            w[i] = 1.0;
            x[i] = 1.0;
            y[i] = 1.0;
        }
        // test z = w * x + y
        Vmath::SIMD::Vvtvp(n, w.data(), x.data(), y.data(), z.data());

        for (size_t i = 0; i < n; ++i)
        {
            BOOST_CHECK_CLOSE(z[i], 2.0, epsilon);
        }

        // ---------------------------------------------------------------------

        // init
        for (size_t i = 0; i < n; ++i)
        {
            w[i] = 0.0;
            x[i] = 1.5;
            y[i] = 0.5;
        }
        // test z = w * x + y
        Vmath::SIMD::Vvtvp(n, w.data(), x.data(), y.data(), z.data());

        for (size_t i = 0; i < n; ++i)
        {
            BOOST_CHECK_CLOSE(z[i], 0.5, epsilon);
        }

        // ---------------------------------------------------------------------

        // init
        for (size_t i = 0; i < n; ++i)
        {
            w[i] = 1.0;
            x[i] = 0.5;
            y[i] = 0.0;
        }
        // test z = w * x + y
        Vmath::SIMD::Vvtvp(n, w.data(), x.data(), y.data(), z.data());

        for (size_t i = 0; i < n; ++i)
        {
            BOOST_CHECK_CLOSE(z[i], 0.5, epsilon);
        }

    }

    BOOST_AUTO_TEST_CASE(TestVvtvvtp)
    {
        using dataType = double;
        constexpr size_t n = 11;
        alignas(tinysimd::simd<dataType>::alignment)
            std::array<dataType, n> v, w, x, y, z;
        dataType epsilon = std::numeric_limits<dataType>::epsilon();

        // init
        for (size_t i = 0; i < n; ++i)
        {
            v[i] = 1.0;
            w[i] = 1.0;
            x[i] = 1.0;
            y[i] = 1.0;
        }
        // test z = v * w + y * z;
        Vmath::SIMD::Vvtvvtp(n, v.data(), w.data(), x.data(), y.data(), z.data());

        for (size_t i = 0; i < n; ++i)
        {
            BOOST_CHECK_CLOSE(z[i], 2.0, epsilon);
        }

        // ---------------------------------------------------------------------

        // init
        for (size_t i = 0; i < n; ++i)
        {
            v[i] = 1.0;
            w[i] = 0.0;
            x[i] = 0.5;
            y[i] = 1.0;
        }
        // test z = v * w + y * z;
        Vmath::SIMD::Vvtvvtp(n, v.data(), w.data(), x.data(), y.data(), z.data());

        for (size_t i = 0; i < n; ++i)
        {
            BOOST_CHECK_CLOSE(z[i], 0.5, epsilon);
        }

        // ---------------------------------------------------------------------

        // init
        for (size_t i = 0; i < n; ++i)
        {
            v[i] = 0.5;
            w[i] = 1.0;
            x[i] = 0.0;
            y[i] = 1.0;
        }
        // test z = v * w + y * z;
        Vmath::SIMD::Vvtvvtp(n, v.data(), w.data(), x.data(), y.data(), z.data());

        for (size_t i = 0; i < n; ++i)
        {
            BOOST_CHECK_CLOSE(z[i], 0.5, epsilon);
        }

    }
}
}
