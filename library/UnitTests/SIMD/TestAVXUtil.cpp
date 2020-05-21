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

#include <AVXOperators/AVXUtil.hpp>

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
namespace AVXUtilsTests
{

    using vec_t = AVX::VecData<double, AVX::SIMD_WIDTH_SIZE>;

    BOOST_AUTO_TEST_CASE(VecData_CopyToAlignedVector_noPad)
    {
        size_t size{AVX::SIMD_WIDTH_SIZE*4};
        Array<OneD, double> aNekArray(size);
        for (size_t i = 0; i < size; ++i)
        {
            aNekArray[i] = i;
        }
        // copy to aligned
        size_t pad = size % AVX::SIMD_WIDTH_SIZE;
        size_t vecSize = size / AVX::SIMD_WIDTH_SIZE + (pad == 0 ? 0 : 1);

        AlignedVector<vec_t> aSIMDVec(vecSize);
        CopyToAlignedVector(aNekArray, aSIMDVec);

        for (size_t i = 0, cnt = 0; i < vecSize; ++i)
        {
            for (size_t j = 0; j < AVX::SIMD_WIDTH_SIZE; ++j, ++cnt)
            {
                BOOST_CHECK_EQUAL(aSIMDVec[i].m_data[j], cnt);
            }
        }

    }

    BOOST_AUTO_TEST_CASE(VecData_CopyToAlignedVector_Pad)
    {
        size_t size{AVX::SIMD_WIDTH_SIZE*4+1};
        Array<OneD, double> aNekArray(size);
        for (size_t i = 0; i < size; ++i)
        {
            aNekArray[i] = i;
        }
        // copy to aligned
        size_t pad = size % AVX::SIMD_WIDTH_SIZE;
        size_t vecSize = size / AVX::SIMD_WIDTH_SIZE + (pad == 0 ? 0 : 1);

        AlignedVector<vec_t> aSIMDVec(vecSize);
        CopyToAlignedVector(aNekArray, aSIMDVec);

        for (size_t i = 0, cnt = 0; i < vecSize; ++i)
        {
            for (size_t j = 0; j < AVX::SIMD_WIDTH_SIZE; ++j, ++cnt)
            {
                if (cnt < size)
                {
                    BOOST_CHECK_EQUAL(aSIMDVec[i].m_data[j], cnt);
                }
                else
                {
                    BOOST_CHECK_EQUAL(aSIMDVec[i].m_data[j], 0);
                }
            }
        }

    }

    BOOST_AUTO_TEST_CASE(VecData_CopyFromAlignedVector_noPad)
    {
        size_t size{AVX::SIMD_WIDTH_SIZE*4};
        Array<OneD, double> aNekArray(size);

        size_t pad = size % AVX::SIMD_WIDTH_SIZE;
        size_t vecSize = size / AVX::SIMD_WIDTH_SIZE + (pad == 0 ? 0 : 1);
        AlignedVector<vec_t> aSIMDVec(vecSize);

        for (size_t i = 0, cnt = 0; i < vecSize; ++i)
        {
            for (size_t j = 0; j < AVX::SIMD_WIDTH_SIZE; ++j, ++cnt)
            {
                aSIMDVec[i].m_data[j] = cnt;
            }
        }

        CopyFromAlignedVector(aSIMDVec, aNekArray);

        for (size_t i = 0; i < size; ++i)
        {
            BOOST_CHECK_EQUAL(aNekArray[i], i);
        }

    }

    BOOST_AUTO_TEST_CASE(VecData_CopyFromAlignedVector_Pad)
    {
        size_t size{AVX::SIMD_WIDTH_SIZE*4+1};
        Array<OneD, double> aNekArray(size);

        size_t pad = size % AVX::SIMD_WIDTH_SIZE;
        size_t vecSize = size / AVX::SIMD_WIDTH_SIZE + (pad == 0 ? 0 : 1);
        AlignedVector<vec_t> aSIMDVec(vecSize);

        for (size_t i = 0, cnt = 0; i < vecSize; ++i)
        {
            for (size_t j = 0; j < AVX::SIMD_WIDTH_SIZE; ++j, ++cnt)
            {
                aSIMDVec[i].m_data[j] = cnt;
            }
        }

        CopyFromAlignedVector(aSIMDVec, aNekArray);

        for (size_t i = 0; i < size; ++i)
        {
            BOOST_CHECK_EQUAL(aNekArray[i], i);
        }

    }




}
}
}