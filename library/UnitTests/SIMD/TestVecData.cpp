///////////////////////////////////////////////////////////////////////////////
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
/// The above copyright notice and this permission notice shall be included
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

#include <AVXOperators/VecData.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>

namespace Nektar
{
namespace AVX
{
namespace VecDataTests
{

    using vec_t = AVX::VecData<double, AVX::SIMD_WIDTH_SIZE>;


    BOOST_AUTO_TEST_CASE(VecData_add)
    {
        double val1 = -4.0;
        double val2 =  2.0;
        double val3 =  2.0;
        vec_t avec1(val1);
        vec_t avec2(val2);
        vec_t avec3(val3);
        vec_t res = avec1 + avec2 * avec3;
        alignas(AVX::SIMD_WIDTH_BYTES) std::array<double, AVX::SIMD_WIDTH_SIZE>
            ascalararr{};
        res.store(ascalararr.data());

        for (size_t i = 0; i < AVX::SIMD_WIDTH_SIZE; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val1 + val2 + val3);
        }

    }

    BOOST_AUTO_TEST_CASE(VecData_store)
    {
        double val = 4.0;
        vec_t avec(val);
        alignas(AVX::SIMD_WIDTH_BYTES) std::array<double, AVX::SIMD_WIDTH_SIZE> ascalararr{};
        avec.store(ascalararr.data());

        for (size_t i = 0; i < AVX::SIMD_WIDTH_SIZE; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val);
        }
    }


    BOOST_AUTO_TEST_CASE(VecData_sqrt)
    {
        double val = 4.0;
        vec_t avec(val);
        vec_t asqrt = sqrt(avec);
        alignas(AVX::SIMD_WIDTH_BYTES) std::array<double, AVX::SIMD_WIDTH_SIZE>
            ascalararr{};
        asqrt.store(ascalararr.data());

        for (size_t i = 0; i < AVX::SIMD_WIDTH_SIZE; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], std::sqrt(val));
        }
    }


    BOOST_AUTO_TEST_CASE(VecData_abs)
    {
        double val = -4.0;
        vec_t avec(val);
        vec_t aabs = abs(avec);
        alignas(AVX::SIMD_WIDTH_BYTES) std::array<double, AVX::SIMD_WIDTH_SIZE>
            ascalararr{};
        aabs.store(ascalararr.data());

        for (size_t i = 0; i < AVX::SIMD_WIDTH_SIZE; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], std::abs(val));
        }

    }

    BOOST_AUTO_TEST_CASE(VecData_load_unload_to_aligned)
    {
        constexpr size_t nDof{5};
        // no padding in load_interleave deinterleave_store
        constexpr size_t nEle{AVX::SIMD_WIDTH_SIZE * 2};
        constexpr size_t nDofBlock = nDof * AVX::SIMD_WIDTH_SIZE;

        constexpr size_t size{nDof*nEle};
        std::array<double,size> dofScalarArr{};
        for (size_t i = 0; i < size; ++i)
        {
            dofScalarArr[i] = i;
        }

        // number of blocks
        size_t nBlock = nEle / AVX::SIMD_WIDTH_SIZE;

        AlignedVector<vec_t> dofVectorArr(nDof);

        double* dataPtr = dofScalarArr.data();
        // loop over blocks AVX::SIMD_WIDTH_SIZE elements at the time
        for (size_t i = 0; i < nBlock; ++i)
        {
            // load
            vec_t::load_interleave(dataPtr, nDof, dofVectorArr);

            // manipulate each block
            for (size_t j = 0; j < nDof; ++j)
            {
                dofVectorArr[j] = dofVectorArr[j] + vec_t(1.0);
            }

            // store
            vec_t::deinterleave_store(dofVectorArr, nDof, dataPtr);
            dataPtr += nDofBlock;
        }

        // check
        for (size_t i = 0; i < size; ++i)
        {
            BOOST_CHECK_EQUAL(dofScalarArr[i], i + 1.0);
        }
    }


}
}
}